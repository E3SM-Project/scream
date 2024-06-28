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


def _limit_sw_down(sfc_sw_down: xr.DataArray, toa_sw_down: xr.DataArray, is_sw_down_toa: xr.DataArray) -> xr.DataArray:
    sfc_sw_down = sfc_sw_down.where(is_sw_down_toa, 0.0)
    # max downwelling based on checking training data
    sfc_sw_down_constraint = toa_sw_down * 0.92  
    sfc_sw_down = sfc_sw_down.where(sfc_sw_down <= sfc_sw_down_constraint, sfc_sw_down_constraint)

    return sfc_sw_down


def _limit_sfc_vis_frac(sfc_vis_frac: xr.DataArray, is_sw_down_toa: xr.DataArray) -> xr.DataArray:
    sfc_vis_frac = sfc_vis_frac.where(is_sw_down_toa, 0.0)
    # sfc vis frac doesn't go lower than ~0.35 where downwelling shortwave exists
    # TODO: should I also add a cap?  Perhaps just add to clip config for training
    vis_frac_ok = cp.logical_and(is_sw_down_toa.data, sfc_vis_frac.data >= 0.35)
    vis_frac_ok = cp.logical_or(vis_frac_ok, cp.logical_not(is_sw_down_toa.data))
    sfc_vis_frac = sfc_vis_frac.where(vis_frac_ok, 0.35)

    return sfc_vis_frac


def _get_constrained_diffusive_fluxes(
    sfc_nir_dif_frac: xr.DataArray,
    sfc_vis_dif_frac: xr.DataArray,
    is_sw_down_toa: xr.DataArray,
) -> xr.DataArray:
    sfc_nir_dif_frac = sfc_nir_dif_frac.where(is_sw_down_toa, 0.0)
    sfc_vis_dif_frac = sfc_vis_dif_frac.where(is_sw_down_toa, 0.0)

    # from notebook, nir dif frac should never be greater than sfc_vis_dif_frac
    nir_vis_dif_inconsistency = sfc_nir_dif_frac > sfc_vis_dif_frac
    sfc_nir_dif_frac = xr.where(
        nir_vis_dif_inconsistency, sfc_vis_dif_frac, sfc_nir_dif_frac
    )
    num_inconsistent_points = nir_vis_dif_inconsistency.sum()
    if num_inconsistent_points.data > 0:
        avg_difference = (sfc_nir_dif_frac - sfc_vis_dif_frac).where(nir_vis_dif_inconsistency).mean()
        logger.info(
            "Number of inconsistent downwelling nir diffuse fractions (> vis diffuse):"
            f" {num_inconsistent_points.data}"
        )
        logger.info(
            "Average difference between downwelling nir and vis diffuse fractions"
            f" at inconsistent points: {avg_difference.data}"
        )
    return sfc_nir_dif_frac, sfc_vis_dif_frac


def get_ML_correction_sfc_fluxes(
    model,
    T_mid,
    qv,
    cos_zenith,
    lat,
    phis,
    sw_flux_dn,
    sfc_alb_dir_vis,
    sfc_alb_dif_vis,
    sfc_alb_dir_nir,
    sfc_alb_dif_nir,
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
    SW_flux_dn_at_model_top = xr.DataArray(sw_flux_dn[:, 0], dims=["ncol"])
    sfc_alb_dir_vis = xr.DataArray(sfc_alb_dir_vis, dims=["ncol"])
    sfc_alb_dif_vis = xr.DataArray(sfc_alb_dif_vis, dims=["ncol"])
    sfc_alb_dir_nir = xr.DataArray(sfc_alb_dir_nir, dims=["ncol"])
    sfc_alb_dif_nir = xr.DataArray(sfc_alb_dif_nir, dims=["ncol"])
    
    ds = xr.Dataset(
        data_vars=dict(
            T_mid=(["ncol", "z"], T_mid),
            qv=(["ncol", "z"], qv),
            lat=(["ncol"], lat),
            surface_geopotential=(["ncol"], phis),
            cos_zenith_angle=(["ncol"], cos_zenith),
        )
    )
    output = predict(model, ds)
    is_sw_down_toa = SW_flux_dn_at_model_top > 0

    # for transmissivity calculation
    atmos_transmissivity = output["shortwave_transmissivity_of_atmospheric_column"].where(is_sw_down_toa, 0.0)
    sfc_sw_down = SW_flux_dn_at_model_top * atmos_transmissivity
    
    # for direct sw_down prediction if used
    # sfc_sw_down = _limit_sw_down(
    #     output["total_sky_downward_shortwave_flux_at_surface"],
    #     SW_flux_dn_at_model_top,
    #     is_sw_down_toa,
    # )

    sfc_vis_frac = _limit_sfc_vis_frac(
        output["downward_vis_fraction_at_surface"],
        is_sw_down_toa,
    )

    sfc_nir_frac = (1 - sfc_vis_frac).where(is_sw_down_toa, 0.0)
    sfc_nir_dif_frac, sfc_vis_dif_frac = _get_constrained_diffusive_fluxes(
        output["downward_nir_diffuse_fraction_at_surface"],
        output["downward_vis_diffuse_fraction_at_surface"],
        is_sw_down_toa,
    )

    sfc_nir_dir_frac = (1 - sfc_nir_dif_frac).where(is_sw_down_toa, 0.0)
    sfc_vis_dir_frac = (1 - sfc_vis_dif_frac).where(is_sw_down_toa, 0.0)

    sfc_flux_vis_dir = sfc_sw_down * sfc_vis_frac * sfc_vis_dir_frac
    sfc_flux_vis_dif = sfc_sw_down * sfc_vis_frac * sfc_vis_dif_frac 
    sfc_flux_nir_dir = sfc_sw_down * sfc_nir_frac * sfc_nir_dir_frac
    sfc_flux_nir_dif = sfc_sw_down * sfc_nir_frac * sfc_nir_dif_frac

    sw_net = (
        sfc_flux_vis_dir * (1 - sfc_alb_dir_vis) +
        sfc_flux_vis_dif * (1 - sfc_alb_dif_vis) +
        sfc_flux_nir_dir * (1 - sfc_alb_dir_nir) +
        sfc_flux_nir_dif * (1 - sfc_alb_dif_nir)
    )

    if (is_sw_down_toa.sum().data > 0):
        logger.info(
            "Min max radiation outputs\n"
            f"\sfc_sw_down_min: {sfc_sw_down.where(is_sw_down_toa).min().data:1.2f}, sfc_sw_down_max: {sfc_sw_down.where(is_sw_down_toa).max().data:1.2f}\n"
            f"\tsw_net_min: {sw_net.where(is_sw_down_toa).min().data:1.2f}, sw_net_max: {sw_net.where(is_sw_down_toa).max().data:1.2f}\n"
            f"\tsfc_vis_frac_min: {sfc_vis_frac.where(is_sw_down_toa).min().data:1.2f}, sfc_vis_frac_max: {sfc_vis_frac.where(is_sw_down_toa).max().data:1.2f}\n"
            f"\tsfc_nir_dif_frac_min: {sfc_nir_dif_frac.where(is_sw_down_toa).min().data:1.2f}, sfc_nir_dif_frac_max: {sfc_nir_dif_frac.where(is_sw_down_toa).max().data:1.2f}\n"
            f"\tsfc_vis_dif_frac_min: {sfc_vis_dif_frac.where(is_sw_down_toa).min().data:1.2f}, sfc_vis_dif_frac_max: {sfc_vis_dif_frac.where(is_sw_down_toa).max().data:1.2f}\n"
        )

    output["sfc_flux_dif_nir"] = sfc_flux_nir_dif
    output["sfc_flux_dif_vis"] = sfc_flux_vis_dif
    output["sfc_flux_dir_nir"] = sfc_flux_nir_dir
    output["sfc_flux_dir_vis"] = sfc_flux_vis_dir
    output["sfc_flux_sw_net"] = sw_net
    output["sfc_flux_sw_dn"] = sfc_sw_down

    return ensure_correction_ordering(output)




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
    sw_flux_dn,
    sfc_alb_dir_vis,
    sfc_alb_dif_vis,
    sfc_alb_dir_nir,
    sfc_alb_dif_nir,
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
    sw_flux_dn = _get_cupy_from_gpu_ptr(sw_flux_dn, (Ncol, Nlev+1), field_dtype_char)
    sfc_alb_dir_vis = _get_cupy_from_gpu_ptr(sfc_alb_dir_vis, (Ncol,), field_dtype_char)
    sfc_alb_dif_vis = _get_cupy_from_gpu_ptr(sfc_alb_dif_vis, (Ncol,), field_dtype_char)
    sfc_alb_dir_nir = _get_cupy_from_gpu_ptr(sfc_alb_dir_nir, (Ncol,), field_dtype_char)
    sfc_alb_dif_nir = _get_cupy_from_gpu_ptr(sfc_alb_dif_nir, (Ncol,), field_dtype_char)
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
            sw_flux_dn,
            sfc_alb_dir_vis,
            sfc_alb_dif_vis,
            sfc_alb_dir_nir,
            sfc_alb_dif_nir,
        )
        sfc_flux_dir_nir[:] = correction_sfc_fluxes["sfc_flux_dir_nir"].data
        sfc_flux_dir_vis[:] = correction_sfc_fluxes["sfc_flux_dir_vis"].data
        sfc_flux_dif_nir[:] = correction_sfc_fluxes["sfc_flux_dif_nir"].data
        sfc_flux_dif_vis[:] = correction_sfc_fluxes["sfc_flux_dif_vis"].data
        sfc_flux_sw_net[:] = correction_sfc_fluxes["sfc_flux_sw_net"].data
        sfc_flux_lw_dn[:] = correction_sfc_fluxes["sfc_flux_lw_dn"].data

    end = time.time()
    logger.debug(
        "ML correction python timing (seconds):  {total: "
        f"{end-start}"
        "}")