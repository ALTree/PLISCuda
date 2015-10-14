#include <cuda_runtime.h>

#include "../../include/cuda/rates.cuh"
#include "../../include/cuda/cuda_utils.cuh"
#include "../../include/cpp/cpp_utils.hpp"

void run_rates_tests();
__global__ void test_react_rates(float * result);
__global__ void test_diff_rates(float * result);
__global__ void test_update_rate_matrix(float * rate_matrix, float * react_rates_array, float * diff_rates_array);
