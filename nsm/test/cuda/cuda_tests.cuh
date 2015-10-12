#include <cuda_runtime.h>

#include "../../include/nsm.cuh"

void run_all();
__global__ void test_react_rates();
__global__ void test_diff_rates();
