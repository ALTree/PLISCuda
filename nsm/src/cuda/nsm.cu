#include "../../include/cuda/nsm.cuh"

__global__ void fill_tau_array(float * tau, int sbc)
{
	int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= sbc)
		return;

	curandState s;
	curand_init(sbi, 0, 0, &s);

	float x = curand_uniform(&s);
	tau[sbi] = x;
	printf("tau[%d] = %f\n", sbi, x);
}

void h_fill_tau_array(thrust::device_vector<float> tau)
{
	fill_tau_array<<<1, tau.size()>>>(thrust::raw_pointer_cast(tau.data()), tau.size());
}
