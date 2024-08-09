import xarray as xr

try:
    import cupy as cp
    import cupy_xarray
    xp = cp
except ImportError:
    import numpy as np
    xp = np

import ml_correction as ml

def test_limit_sw_down():
    # should be no shortwave at surface where there is no shortwave at TOA
    # sfc shortwave shouldn't be more than maximium of ~0.92*TOA shortwave
    toa_sw_down = xr.DataArray(xp.array([0, 10, 10], dtype=xp.float32), dims=["ncol"])
    sw_down_mask = xr.DataArray(xp.array([0, 1, 1], dtype=xp.bool_), dims=["ncol"])
    sfc_sw_down = xr.DataArray(xp.array([5, 11, 3], dtype=xp.float32), dims=["ncol"])

    constrained = ml._limit_sw_down(sfc_sw_down, toa_sw_down, sw_down_mask)
    expected = xp.array([0, 9.2, 3], dtype=xp.float32)
    xp.testing.assert_array_equal(constrained.data, expected)


def test_limit_sfc_vis_frac():
    # shouldn't be lower than 0.35
    # should be 0 where there is no shortwave down at TOA
    sw_down_mask = xr.DataArray(xp.array([0, 1, 1], dtype=xp.bool_), dims=["ncol"])
    sfc_vis_frac = xr.DataArray(xp.array([0.3, 0.5, 0.2], dtype=xp.float32), dims=["ncol"])

    constrained = ml._limit_sfc_vis_frac(sfc_vis_frac, sw_down_mask)
    expected = xp.array([0.0, 0.5, 0.35], dtype=xp.float32)
    xp.testing.assert_array_almost_equal(constrained.data, expected)


def test_get_constrained_diffusive_fluxes():
    # the nir diffusive should always be less than or equal to the vis diffusive
    vis_diffusive = xr.DataArray(xp.array([0.5, 0.2], dtype=xp.float32), dims=["ncol"])
    nir_diffusive = xr.DataArray(xp.array([0.2, 0.6], dtype=xp.float32), dims=["ncol"])

    constrained_nir_dif, constrained_vis_dif = ml._get_constrained_diffusive_fluxes(
        nir_diffusive, vis_diffusive
    )
    expected_nir_dif = xp.array([0.2, 0.2], dtype=xp.float32)
    xp.testing.assert_array_equal(constrained_nir_dif.data, expected_nir_dif)
    expected_vis_dif = xp.array([0.5, 0.2], dtype=xp.float32)
    xp.testing.assert_array_equal(constrained_vis_dif.data, expected_vis_dif)


def test_get_direct_diffusive_fractions():
    # should be 0 where there is no shortwave down at TOA
    # should default to direct when both are provided in the dataset
    # check that providing either or also works
    # check that NIR diffusive is constrained by VIS diffusive

    is_sw_down_toa = xr.DataArray(xp.array([0, 1, 1], dtype=bool), dims=["ncol"])
    diffuse = xr.DataArray(xp.array([1.0, 0.4, 0.7], dtype=xp.float32), dims=["ncol"])
    direct = xr.DataArray(xp.array([1.0, 0.6, 0.3], dtype=xp.float32), dims=["ncol"])

    direct_ds = xr.Dataset({
        ml.VIS_DIRECT_KEY: direct,
        ml.NIR_DIRECT_KEY: direct,
    })

    diffusive_ds = xr.Dataset({
        ml.NIR_DIFFUSE_KEY: diffuse,
        ml.VIS_DIFFUSE_KEY: diffuse,
    })
    for ds in [direct_ds, diffusive_ds]:
        vis_dir, vis_dif, nir_dir, nir_dif = ml._get_direct_diffusive_fractions(ds, is_sw_down_toa)
        xp.testing.assert_array_almost_equal(vis_dir.data, (direct * is_sw_down_toa).data)
        xp.testing.assert_array_almost_equal(vis_dif.data, (diffuse * is_sw_down_toa).data)
        xp.testing.assert_array_almost_equal(nir_dir.data, (direct * is_sw_down_toa).data)
        xp.testing.assert_array_almost_equal(nir_dif.data, (diffuse * is_sw_down_toa).data)


def test_get_direct_diffusive_fractions_all_inputs_provided():
    # should be 0 where there is no shortwave down at TOA
    # should default to direct when both are provided in the dataset
    # check that providing either or also works
    # check that NIR diffusive is constrained by VIS diffusive

    is_sw_down_toa = xr.DataArray(xp.array([0, 1, 1], dtype=bool), dims=["ncol"])
    diffuse = xr.DataArray(xp.array([1.0, 0.5, 0.2], dtype=xp.float32), dims=["ncol"])
    direct = xr.DataArray(xp.array([1.0, 0.6, 0.3], dtype=xp.float32), dims=["ncol"])
    expected_diffuse = xr.DataArray(xp.array([0.0, 0.4, 0.7], dtype=xp.float32), dims=["ncol"])

    ds = xr.Dataset({
        ml.VIS_DIRECT_KEY: direct,
        ml.NIR_DIRECT_KEY: direct,
        ml.NIR_DIFFUSE_KEY: diffuse,
        ml.VIS_DIFFUSE_KEY: diffuse
    })

    vis_dir, vis_dif, nir_dir, nir_dif = ml._get_direct_diffusive_fractions(ds, is_sw_down_toa)
    xp.testing.assert_array_almost_equal(vis_dir.data, (direct * is_sw_down_toa).data)
    xp.testing.assert_array_almost_equal(vis_dif.data, (expected_diffuse).data)
    xp.testing.assert_array_almost_equal(nir_dir.data, (direct * is_sw_down_toa).data)
    xp.testing.assert_array_almost_equal(nir_dif.data, (expected_diffuse).data)


def test_get_direct_diffusive_fractions_diffuse_constraint():
    # check that diffusive constraint for NIR is applied and VIS is updated

    is_sw_down_toa = xr.DataArray(xp.array([0, 1, 1], dtype=bool), dims=["ncol"])
    diffuse_nir = xr.DataArray(xp.array([1.0, 0.5, 0.2], dtype=xp.float32), dims=["ncol"])
    diffuse_vis = xr.DataArray(xp.array([1.0, 0.2, 0.2], dtype=xp.float32), dims=["ncol"])
    expected_diffuse = xr.DataArray(xp.array([0.0, 0.2, 0.2], dtype=xp.float32), dims=["ncol"])
    expected_direct = xr.DataArray(xp.array([0.0, 0.8, 0.8], dtype=xp.float32), dims=["ncol"])

    ds = xr.Dataset({
        ml.NIR_DIFFUSE_KEY: diffuse_nir,
        ml.VIS_DIFFUSE_KEY: diffuse_vis
    })

    vis_dir, vis_dif, nir_dir, nir_dif = ml._get_direct_diffusive_fractions(ds, is_sw_down_toa)
    xp.testing.assert_array_almost_equal(vis_dir.data, expected_direct.data)
    xp.testing.assert_array_almost_equal(vis_dif.data, expected_diffuse.data)
    xp.testing.assert_array_almost_equal(nir_dir.data, expected_direct.data)
    xp.testing.assert_array_almost_equal(nir_dif.data, expected_diffuse.data)

    