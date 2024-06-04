import numpy as np
import cupy as cp
import xarray as xr
import cupy_xarray
import datetime
import logging
import time
from vcm import cos_zenith_angle
from scream_run.steppers.machine_learning import (
    MachineLearningConfig,
    open_model,
    predict,
    predict_with_qv_constraint
)

logger = logging.getLogger(__name__)

TQ = "tq"
UV = "uv"
FLUX = "flux"

ML_MODELS = {
    TQ: None,
    UV: None,
    FLUX: None,
}


def _load_model(key, path):
    # No need to pass model objects back and forth into C++
    # Just store the model in a global dict
    global ML_MODELS
    model = ML_MODELS[key]
    
    if model is None and path.lower() != "none":
        logger.info(f"Loading ML {key} model from {path}")
        config = MachineLearningConfig(models=[path])
        model = open_model(config)
        ML_MODELS[key] = model
    elif path.lower() == "none":
        logger.debug("skipping ML {key} model, no path provided")
    else:
        logger.debug(f"skipping ML {key} model, already loaded")
    return model


def load_all_models(tq_model_path, uv_model_path, flx_model_path):
    """Load all models into the global model dict"""
    _load_model(TQ, tq_model_path)
    _load_model(UV, uv_model_path)
    _load_model(FLUX, flx_model_path)


def initialize_logging(level=logging.INFO):
    logging.basicConfig(level=level)


def ensure_correction_ordering(correction):
    """Ensure that the ordering of the correction is always (ncol, z)"""
    for key in correction:
        if "z" in correction[key].dims:
            correction[key] = correction[key].transpose("ncol", "z")
    return correction


def get_ML_correction_dQ1_dQ2(model, T_mid, qv, cos_zenith, lat, phis, dt):
    """Get ML correction for air temperature (dQ1) and specific humidity (dQ2)

    Args:
        model: pre-trained ML model for dQ1 and dQ2
        T_mid: air temperature
        qv: specific humidity
        cos_zenith: cosine zenith angle
        dt: time step (s)
    """
    # TODO: standardize reference DataArray dims 1D or 2D
    ds = xr.Dataset(
        data_vars=dict(
            T_mid=(["ncol", "z"], T_mid),
            qv=(["ncol", "z"], qv),
            cos_zenith_angle=(["ncol"], cos_zenith),
            lat=(["ncol"], lat),
            surface_geopotential=(["ncol"], phis),
        )
    )
    return ensure_correction_ordering(predict_with_qv_constraint(model, ds, dt))


def get_ML_correction_dQu_dQv(model, T_mid, qv, cos_zenith, lat, phis, u, v, dt):
    """Get ML correction for eastward wind (dQu or dQxwind) and northward wind (dQv or dQywind)

    Args:
        model: pre-trained ML model for dQu and dQv
        T_mid: air  temperature
        qv: specific humidity
        cos_zenith: cosine zenith angle
        lat: latitude
        phis: surface geopotential
        u: horizontal wind in x-direction
        v: horizontal wind in y-direction
        dt: time step (s)
    """
    ds = xr.Dataset(
        data_vars=dict(
            T_mid=(["ncol", "z"], T_mid),
            qv=(["ncol", "z"], qv),
            U=(["ncol", "z"], u),
            V=(["ncol", "z"], v),
            lat=(["ncol"], lat),
            surface_geopotential=(["ncol"], phis),
            cos_zenith_angle=(["ncol"], cos_zenith),
        )
    )
    output = ensure_correction_ordering(predict(model, ds))
    # rename dQxwind and dQywind to dQu and dQv if needed
    if "dQxwind" in output.keys():
        output["dQu"] = output.pop("dQxwind")
    if "dQywind" in output.keys():
        output["dQv"] = output.pop("dQywind")

    return output


def get_ML_correction_sfc_fluxes(
    model,
    T_mid,
    qv,
    cos_zenith,
    lat,
    phis,
    # sfc_alb_dif_vis,
    # sw_flux_dn,
):    
    """Get ML correction for overriding surface fluxes (net shortwave and downward longwave)
    ML model should have the following output variables:
        net_shortwave_sfc_flux_via_transmissivity
        override_for_time_adjusted_total_sky_downward_longwave_flux_at_surface

    Args:
        model: pre-trained ML model for radiative fluxes
        T_mid: air  temperature
        qv: specific humidity
        cos_zenith: cosine zenith angle
        lat: latitude
        phis: surface geopotential
        sfc_alb_dif_vis: surface albedo for diffuse shortwave radiation
        sw_flux_dn: downward shortwave flux
        dt: time step (s)
    """    
    # SW_flux_dn_at_model_top = sw_flux_dn[:, 0]
    ds = xr.Dataset(
        data_vars=dict(
            T_mid=(["ncol", "z"], T_mid),
            qv=(["ncol", "z"], qv),
            lat=(["ncol"], lat),
            surface_geopotential=(["ncol"], phis),
            cos_zenith_angle=(["ncol"], cos_zenith),
            # surface_diffused_shortwave_albedo=(["ncol"], sfc_alb_dif_vis),
            # total_sky_downward_shortwave_flux_at_top_of_atmosphere=(
            #     ["ncol"],
            #     SW_flux_dn_at_model_top,
            # ),
        )
    )
    return predict(model, ds)


def _get_cupy_from_gpu_ptr(ptr, shape, dtype_char):
    dtype = cp.dtype(dtype_char)
    total_size = np.product(shape) * dtype.itemsize
    data_from_ptr = cp.ndarray(
        shape, 
        dtype=dtype,
        memptr=cp.cuda.MemoryPointer(
            cp.cuda.UnownedMemory(ptr, total_size, None), 0
        )
    )
    return data_from_ptr


def update_fields(
    field_dtype_char,
    qv,
    T_mid,
    u,
    v,
    lat,
    lon,
    phis,
    sfc_flux_dir_nir,
    sfc_flux_dir_vis,
    sfc_flux_dif_nir,
    sfc_flux_dif_vis,
    sfc_flux_sw_net,
    sfc_flux_lw_dn,
    Ncol,
    Nlev,
    num_tracers,
    dt,
    current_time,
    tq_model_path,
    uv_model_path,
    flux_model_path,
):
    """
    T_mid: temperature
    qv: specific humidity
    u: x-component of wind
    v: y-component of wind
    lat: latitude
    lon: longitude
    phis: surface geopotential
    SW_flux_dn_at_model_top: downwelling shortwave flux at the top of the model
    sfc_alb_dif_vis: surface diffuse shortwave albedo
    sfc_flux_sw_net
    sfc_flux_lw_dn
    Ncol: number of columns
    Nlev: number of levels
    num_tracers: number of tracers
    dt: time step (s)
    current_time: current time in the format "YYYY-MM-DD HH:MM:SS"
    tq_model_path: path to the temperature and specific humidity ML model
    uv_model_path: path to the wind ML model
    flux_model_path: path to the surface fluxes ML model
    """
    start = time.time()
    qv = _get_cupy_from_gpu_ptr(qv, (Ncol, num_tracers, Nlev), field_dtype_char)
    T_mid = _get_cupy_from_gpu_ptr(T_mid, (Ncol, Nlev), field_dtype_char)
    u = _get_cupy_from_gpu_ptr(u, (Ncol, Nlev), field_dtype_char)
    v = _get_cupy_from_gpu_ptr(v, (Ncol, Nlev), field_dtype_char)
    lat = _get_cupy_from_gpu_ptr(lat, (Ncol,), field_dtype_char)
    lon = _get_cupy_from_gpu_ptr(lon, (Ncol,), field_dtype_char)
    phis = _get_cupy_from_gpu_ptr(phis, (Ncol,), field_dtype_char)
    # sw_flux_dn = _get_cupy_from_gpu_ptr(sw_flux_dn, (Ncol, Nlev+1), field_dtype_char)
    # sfc_alb_dif_vis = _get_cupy_from_gpu_ptr(sfc_alb_dif_vis, (Ncol,), field_dtype_char)
    sfc_flux_dir_nir = _get_cupy_from_gpu_ptr(sfc_flux_dir_nir, (Ncol,), field_dtype_char)
    sfc_flux_dir_vis = _get_cupy_from_gpu_ptr(sfc_flux_dir_vis, (Ncol,), field_dtype_char)
    sfc_flux_dif_nir = _get_cupy_from_gpu_ptr(sfc_flux_dif_nir, (Ncol,), field_dtype_char)
    sfc_flux_dif_vis = _get_cupy_from_gpu_ptr(sfc_flux_dif_vis, (Ncol,), field_dtype_char)
    sfc_flux_sw_net = _get_cupy_from_gpu_ptr(sfc_flux_sw_net, (Ncol,), field_dtype_char)
    sfc_flux_lw_dn = _get_cupy_from_gpu_ptr(sfc_flux_lw_dn, (Ncol,), field_dtype_char)
    current_datetime = datetime.datetime.strptime(current_time, "%Y-%m-%d %H:%M:%S")
    cos_zenith = cos_zenith_angle(
        current_datetime,
        lon,
        lat,
    )

    # Kokkos crashes if I import tensorflow too early 
    # (i.e., during MLCorrection::initialize_impl)
    # so defer the import until it's needed for models here.
    # models are loaded into global state so they are only loaded once
    # despite being called each time we ask for a correction
    model_tq = _load_model(TQ, tq_model_path)
    model_uv = _load_model(UV, uv_model_path)
    model_sfc_fluxes = _load_model(FLUX, flux_model_path)
    
    device = qv.device.id
    with cp.cuda.Device(device):
        qv_0 = cp.ascontiguousarray(qv[:, 0, :])
    
    if model_tq is not None:
        correction_tq = get_ML_correction_dQ1_dQ2(
            model_tq, 
            T_mid, 
            qv_0, 
            cos_zenith,
            lat,
            phis,
            dt
        )
        T_mid[:, :] += correction_tq["dQ1"].data * dt
        qv[:, 0, :] += correction_tq["dQ2"].data * dt
    if model_uv is not None:
        correction_uv = get_ML_correction_dQu_dQv(
            model_uv, T_mid, qv_0, cos_zenith, lat, phis, u, v, dt
        )
        u[:, :] += correction_uv["dQu"].data * dt
        v[:, :] += correction_uv["dQv"].data * dt
    if model_sfc_fluxes is not None:
        correction_sfc_fluxes = get_ML_correction_sfc_fluxes(
            model_sfc_fluxes,
            T_mid,
            qv_0,
            cos_zenith,
            lat,
            phis,
            # sfc_alb_dif_vis,
            # sw_flux_dn,
        )
        # sfc_flux_sw_net[:] = correction_sfc_fluxes["net_shortwave_sfc_flux_via_transmissivity"].data
        # sfc_flux_lw_dn[:] = correction_sfc_fluxes["override_for_time_adjusted_total_sky_downward_longwave_flux_at_surface"].data
        sfc_flux_dir_nir[:] = correction_sfc_fluxes["sfc_flux_dir_nir"].data
        sfc_flux_dir_vis[:] = correction_sfc_fluxes["sfc_flux_dir_vis"].data
        sfc_flux_dif_nir[:] = correction_sfc_fluxes["sfc_flux_dif_nir"].data
        sfc_flux_dif_vis[:] = correction_sfc_fluxes["sfc_flux_dif_vis"].data
        sfc_flux_sw_net[:] = correction_sfc_fluxes["sfc_flux_sw_net"].data
        sfc_flux_lw_dn[:] = correction_sfc_fluxes["sfc_flux_lw_dn"].data

    end = time.time()
    logger.info(
        "ML correction python timing (seconds):  {total: "
        f"{end-start}"
        "}")