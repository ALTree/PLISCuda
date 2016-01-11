#ifndef LOG_CUH_
#define LOG_CUH_

#include <cuda_runtime.h>

#include "cuda_utils.cuh"
#include "constants.cuh"

__global__ void log_data(int * state, bool * sbv_to_log, bool * spc_to_log, int * log_data);

#endif /* LOG_CUH_ */
