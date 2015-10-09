#include <cuda_runtime.h>

#include "../../include/nsm.cuh"

void run_all();

__global__ void driver();
__device__ void test_react_rates();
