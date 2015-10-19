#ifndef NSM_CUH_
#define NSM_CUH_

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>


__global__ void fill_tau_array(float * tau, int sbc);

void h_fill_tau_array(thrust::device_vector<float> tau);


#endif /* NSM_CUH_ */
