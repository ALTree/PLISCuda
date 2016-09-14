#ifndef LOG_CUH_
#define LOG_CUH_

#include <cuda_runtime.h>

#include "constants.cuh"

__global__ void log_data(int * state, unsigned int * sbv_to_log, int sbv_to_log_len, bool * spc_to_log,
						 int spc_to_log_len, int * log_data);

#endif /* LOG_CUH_ */