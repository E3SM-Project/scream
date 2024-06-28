import xarray as xr
import cupy as cp
import cupy_xarray

from ml_correction import _limit_sw_down, _limit_sfc_vis_frac, _get_constrained_diffusive_fluxes

def test_limit_sw_down():
    # should be no shortwave at surface where there is no shortwave at TOA
    # sfc shortwave shouldn't be more than maximium of ~0.92*TOA shortwave
    toa_sw_down = xr.DataArray(cp.array([0, 10, 10], dtype=cp.float32), dims=["ncol"])
    sw_down_mask = xr.DataArray(cp.array([0, 1, 1], dtype=cp.bool_), dims=["ncol"])
    sfc_sw_down = xr.DataArray(cp.array([5, 11, 3], dtype=cp.float32), dims=["ncol"])

    constrained = _limit_sw_down(sfc_sw_down, toa_sw_down, sw_down_mask)
    expected = cp.array([0, 9.2, 3], dtype=cp.float32)
    cp.testing.assert_array_equal(constrained.data, expected)


def test_limit_sfc_vis_frac():
    # shouldn't be lower than 0.35
    # should be 0 where there is no shortwave down at TOA
    sw_down_mask = xr.DataArray(cp.array([0, 1, 1], dtype=cp.bool_), dims=["ncol"])
    sfc_vis_frac = xr.DataArray(cp.array([0.3, 0.5, 0.2], dtype=cp.float32), dims=["ncol"])

    constrained = _limit_sfc_vis_frac(sfc_vis_frac, sw_down_mask)
    expected = cp.array([0.0, 0.5, 0.35], dtype=cp.float32)
    cp.testing.assert_array_almost_equal(constrained.data, expected)


def test_get_constrained_diffusive_fluxes():
    # should be 0 where there is no shortwave down at TOA
    # the nir diffusive should always be less than or equal to the vis diffusive
    sw_down_mask = xr.DataArray(cp.array([0, 1, 1], dtype=cp.bool_), dims=["ncol"])
    vis_diffusive = xr.DataArray(cp.array([0.3, 0.5, 0.2], dtype=cp.float32), dims=["ncol"])
    nir_diffusive = xr.DataArray(cp.array([0.4, 0.2, 0.6], dtype=cp.float32), dims=["ncol"])

    constrained_nir_dif, constrained_vis_dif = _get_constrained_diffusive_fluxes(
        nir_diffusive, vis_diffusive, sw_down_mask
    )
    expected_nir_dif = cp.array([0.0, 0.2, 0.2], dtype=cp.float32)
    cp.testing.assert_array_equal(constrained_nir_dif.data, expected_nir_dif)
    expected_vis_dif = cp.array([0.0, 0.5, 0.2], dtype=cp.float32)
    cp.testing.assert_array_equal(constrained_vis_dif.data, expected_vis_dif)
    