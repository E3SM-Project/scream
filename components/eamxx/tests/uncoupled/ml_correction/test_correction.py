import numpy as np
import cupy as cp
import h5py
import xarray as xr
import cupy_xarray
import fv3fit
from scream_run.steppers.machine_learning import (
    MultiModelAdapter,
    predict,
)


def get_ML_model(model_path):
    if model_path == "NONE":
        return None
    config = MachineLearningConfig(models=[model_path])
    model = open_model(config)
    return model


def sample_ML_prediction(
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
    output = predict(model, input_data, dt=1.0)
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
    prediction = sample_ML_prediction(Nlev, data_from_ptr[1, :], model_tq, model_uv)
    data_from_ptr[1, :] = prediction


def modify_view(data, Ncol, Nlev, model_tq, model_uv):
    # data comes in as 1D Numpy array, a view constructed from C++ memory
    data_cp = cp.asarray(data)
    data_np = cp.asnumpy(data_cp)
    np.testing.assert_array_equal(data_np, data)
    data = data.reshape((-1, Nlev))
    data_cp_rshp = cp.reshape(data_cp, (-1, Nlev))
    prediction = sample_ML_prediction(Nlev, data_cp_rshp[1, :], model_tq, model_uv)
    data[1, :] = prediction.get()


if __name__ == "__main__":
    Ncol = 2
    Nlev = 72
    model_tq = "NONE"
    model_uv = "NONE"
    data = cp.random.rand(Ncol, Nlev)
    modify_view(data, Ncol, Nlev, model_tq, model_uv)
    print(data)
