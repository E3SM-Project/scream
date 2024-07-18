import numpy as np

try:
    import cupy as cp
    import cupy_xarray
except ImportError:
    cp = np
import h5py
import xarray as xr
import fv3fit
import logging
from scream_run.steppers.machine_learning import (
    MachineLearningConfig,
    MultiModelAdapter,
    predict,
    open_model,
)
from typing import Union

logger = logging.getLogger(__name__)


def get_ML_model(model_path):
    if model_path == "NONE":
        return None
    config = MachineLearningConfig(models=[model_path])
    model = open_model(config)
    return model


def sample_ML_prediction(nz, input_data, ML_model_tq: fv3fit.PureKerasModel):
    """
    This function is used to generate a sample ML prediction for the given input data.
    We use a constant output predictor to generate the prediction.
    """
    if isinstance(input_data, cp.ndarray):
        inputs = {
            "T_mid": (["ncol", "z"], input_data),
            "qv": (["ncol", "z"], input_data),
            "cos_zenith_angle": (["ncol"], cp.full((input_data.shape[0]), 0.5)),
        }
    else:
        inputs = {
            "T_mid": (["ncol", "z"], input_data),
            "qv": (["ncol", "z"], input_data),
            "cos_zenith_angle": (["ncol"], np.full((input_data.shape[0]), 0.5)),
        }

    input_data = xr.Dataset(inputs)
    output = predict(ML_model_tq, input_data)
    # just check the output from one variable to see if it works
    return output["dQ1"].data


def sample_ML_prediction_dummy(
    nz: int,
    input_data,
    ML_model_tq: str,
    ML_model_uv: str,
    ML_uses_GPU: bool,
):
    """
    This function is used to generate a sample ML prediction for the given input data.
    We use a constant output predictor to generate the prediction.
    """
    output_variables = ["qv"]
    if ML_uses_GPU:
        outputs = {
            "qv": cp.full(nz, 1e-4),
        }
    else:
        outputs = {
            "qv": np.full(nz, 1e-4),
        }
    predictor = fv3fit.testing.ConstantOutputPredictor(
        input_variables=["qv"],
        output_variables=output_variables,
    )
    predictor.set_outputs(**outputs)
    model = MultiModelAdapter([predictor])
    if len(input_data.shape) < 2:
        if ML_uses_GPU:
            input_data = input_data[cp.newaxis, :]
        else:
            input_data = input_data[np.newaxis, :]
    input_data = xr.Dataset({"qv": xr.DataArray(data=input_data, dims=["ncol", "z"])})
    output = predict(model, input_data)
    return output["qv"].data


def modify_view_gpu(ptr, dtype_char, Ncol, Nlev, model_tq, model_uv):
    # data comes in as 1D Numpy array, a view constructed from C++ memory
    dtype = cp.dtype(dtype_char)
    data_from_ptr = cp.ndarray(
        (Ncol, Nlev),
        dtype=dtype,
        memptr=cp.cuda.MemoryPointer(
            cp.cuda.UnownedMemory(ptr, Ncol * Nlev * dtype.itemsize, None), 0
        ),
    )
    prediction = sample_ML_prediction_dummy(
        Nlev, data_from_ptr, model_tq, model_uv, True
    )
    data_from_ptr[1, :] = prediction[1, :]


def build_gpu_array(ptr, dtype_char, Ncol, Nlev):
    dtype = cp.dtype(dtype_char)
    data_from_ptr = cp.ndarray(
        (Ncol, Nlev),
        dtype=dtype,
        memptr=cp.cuda.MemoryPointer(
            cp.cuda.UnownedMemory(ptr, Ncol * Nlev * dtype.itemsize, None), 0
        ),
    )
    return data_from_ptr


def compare_cpu_gpu_arrays(cpu_data, gpu_ptr, dtype_char, Ncol, Nlev):
    cpu = cpu_data.reshape((-1, Nlev))
    assert type(cpu) == np.ndarray
    gpu = build_gpu_array(gpu_ptr, dtype_char, Ncol, Nlev)
    assert type(gpu) == cp.ndarray
    np.testing.assert_allclose(cpu, gpu.get())


def modify_view_cpu(data, Ncol, Nlev, model_tq, model_uv):
    # data comes in as 1D Numpy array, a view constructed from C++ memory
    data = data.reshape((-1, Nlev))
    prediction = sample_ML_prediction_dummy(Nlev, data, model_tq, model_uv, False)
    data[1, :] = prediction[1, :]
