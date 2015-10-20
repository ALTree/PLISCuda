#ifndef NSM_CUH_
#define NSM_CUH_

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>


__global__ void fill_tau_array(float * tau, int sbc);

void h_fill_tau_array(thrust::device_vector<float> &tau);
int h_get_min_tau(thrust::device_vector<float> &tau);


#endif /* NSM_CUH_ */
