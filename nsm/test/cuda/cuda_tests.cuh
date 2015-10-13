#include <cuda_runtime.h>

#include "../../include/cuda/rates.cuh"
#include "../../include/cuda/cuda_utils.cuh"
#include "../../include/cpp/cpp_utils.hpp"

void run_all();
__global__ void test_react_rates(float * result);
__global__ void test_diff_rates(float * result);
// __global__ void test_rate_matrix_row(float *);
