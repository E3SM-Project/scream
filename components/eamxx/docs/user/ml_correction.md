# ML Correction in EAMxx

Machine learning correction is used to emulate lower resolution simulation to behave like its higher resolution counterpart. This project is under development and is only supported on selected machines.

## Example setup

To enable ML correction as a process:

```shell
./atmchange physics::atm_procs_list="mac_aero_mic,rrtmgp,mlcorrection"
```

The following options can be specified:

```shell
./atmchange ML_model_path_tq=/path/to/pretrained/temperature_and_specific_humidity_model
./atmchange ML_model_path_temperature=/path/to/pretrained/temperature_only_model
./atmchange ML_model_path_uv=/path/to/pretrained/u_and_v_model
./atmchange ML_model_path_sfc_fluxes=/path/to/pretrained/surface_fluxes_model
```

## Supported Machines

- Ruby: CPU only machine. Shared python environment can be loaded with `source /usr/WS1/climdat/python_venv/3.9.2/screamML/bin/activate`
- Perlmutter: Both CPU and GPU supported. Shared python environment located in `/global/common/software/m4492/fv3net-shared-py39`. Additional modules and flags need to be specified:

```shell
module load python
conda activate /global/common/software/m4492/fv3net-shared-py39
export TF_USE_LEGACY_KERAS=1
# The followings for running on GPUs:
module load cudnn/8.9.3_cuda12
export TF_FORCE_GPU_ALLOW_GROWTH=true
export XLA_FLAGS=--xla_gpu_cuda_data_dir=/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/cuda/12.2

```
