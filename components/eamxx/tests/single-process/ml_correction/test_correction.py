import numpy as np
import cupy as cp
import h5py
import xarray as xr
import cupy_xarray
import tensorflow as tf
import fv3fit
import logging
from scream_run.steppers.machine_learning import (
    MachineLearningConfig,
    MultiModelAdapter,
    predict,
    open_model
)
from typing import Union

logger = logging.getLogger(__name__)


def get_ML_model(model_path):
    if model_path == "NONE":
        return None
    config = MachineLearningConfig(models=[model_path])
    model = open_model(config)
    return model


def sample_ML_prediction(nz, input_data: Union[cp.ndarray, np.ndarray], ML_model_tq: fv3fit.PureKerasModel):
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
    nz: int, input_data: cp.ndarray, ML_model_tq: str, ML_model_uv: str
):
    """
    This function is used to generate a sample ML prediction for the given input data.
    We use a constant output predictor to generate the prediction.
    """
    output_variables = ["qv"]
    outputs = {
        "qv": cp.full(nz, 1e-4),
    }
    predictor = fv3fit.testing.ConstantOutputPredictor(
        input_variables=["qv"],
        output_variables=output_variables,
    )
    predictor.set_outputs(**outputs)
    model = MultiModelAdapter([predictor])
    if len(input_data.shape) < 2:
        input_data = input_data[cp.newaxis, :]
    input_data = xr.Dataset({"qv": xr.DataArray(data=input_data, dims=["ncol", "z"])})
    output = predict(model, input_data)
    return output["qv"].data


def test_ptr(ptr, data, Ncol, Nlev):
    # TODO: figure out how to move data in C++ onto a GPU and send it
    #       and also send original np array to check data against
    print(f"Pointer received: {ptr}")
    print(f"Data received: {type(data)}")
    print(f"Ncol: {Ncol}, Nlev: {Nlev}")
    data_from_ptr = cp.ndarray(
        (Ncol, Nlev), 
        dtype=cp.float64,
        memptr=cp.cuda.MemoryPointer(
            cp.cuda.UnownedMemory(ptr, Ncol * Nlev * 8, None), 0
        )
    )
    print(f"Data from pointer device: {data_from_ptr.device}")
    print(f"Data type from pointer: {data_from_ptr.dtype}")
    print(f"Data type from host (original): {data.dtype}")
    print(f"Data from pointer: {data_from_ptr}")
    print(f"Data from host (original): {data}")
    pass


def modify_view_gpu(ptr, dtype_char, Ncol, Nlev, model_tq, model_uv):
    # data comes in as 1D Numpy array, a view constructed from C++ memory
    dtype = cp.dtype(dtype_char)
    data_from_ptr = cp.ndarray(
        (Ncol, Nlev), 
        dtype=dtype,
        memptr=cp.cuda.MemoryPointer(
            cp.cuda.UnownedMemory(ptr, Ncol * Nlev * dtype.itemsize, None), 0
        )
    )
    prediction = sample_ML_prediction_dummy(Nlev, data_from_ptr, model_tq, model_uv)
    data_from_ptr[1, :] = prediction[1, :]


def modify_view(data, Ncol, Nlev, model_tq, model_uv):
    # data comes in as 1D Numpy array, a view constructed from C++ memory
    data_cp = cp.asarray(data)
    data = data.reshape((-1, Nlev))
    data_cp_rshp = cp.reshape(data_cp, (-1, Nlev))
    prediction = sample_ML_prediction_dummy(Nlev, data_cp_rshp, model_tq, model_uv)
    data[1, :] = prediction[1, :].get()


def gpu_handoff_real_model(ptr, dtype_char, Ncol, Nlev, model_path):
    dtype = cp.dtype(dtype_char)
    data_from_ptr = cp.ndarray(
        (Ncol, Nlev), 
        dtype=dtype,
        memptr=cp.cuda.MemoryPointer(
            cp.cuda.UnownedMemory(ptr, Ncol * Nlev * dtype.itemsize, None), 0
        )
    )
    model = get_ML_model(model_path)
    gpu_result = sample_ML_prediction(Nlev, data_from_ptr, model).get()

    with tf.device("/cpu:0"):
        cpu_data = data_from_ptr.get()
        logger.info(f"CPU data type: {type(cpu_data)}")
        cpu_result = sample_ML_prediction(Nlev, cpu_data, model)

    cp.testing.assert_array_almost_equal(cpu_result, gpu_result, decimal=6)


def gpu_handoff_real_model_local_test(data, Ncol, Nlev, model_path):
    model = get_ML_model(model_path)
    gpu_result = sample_ML_prediction(Nlev, data, model).get()

    with tf.device("/cpu:0"):
        cpu_data = data.get()
        cpu_result = sample_ML_prediction(Nlev, cpu_data, model)

    cp.testing.assert_array_almost_equal(cpu_result, gpu_result, decimal=6)


if __name__ == "__main__":
    Ncol = 2
    Nlev = 128
    model_tq_path = "/global/cfs/cdirs/m4492/corrective_ml/case_ne120_to_ne30_20240422_novertical_remap/2024-04-19-tq-test-1cc4324d6138/no-tapering"
    data = np.random.rand(Ncol, Nlev)
    gpu_handoff_real_model_local_test(cp.asarray(data), Ncol, Nlev, model_tq_path)
