#include "../../include/cuda/nsm.cuh"

__device__ int choose_rand_reaction(int sbc, int rc, float * rate_matrix, float * react_rates_array, float rand)
{
	int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= sbc)
		return -1;

	float sum = rate_matrix[sbi];
	float scaled_sum = sum * rand;
	float partial_sum = 0;

	int ri = 0;
	while(partial_sum <= scaled_sum) {
		partial_sum += react_rates_array[ri * sbc + sbi];
		ri++;
	}

	return ri-1;
}

__global__ void fill_tau_array(float * tau, int sbc)
{
	int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= sbc)
		return;

	curandState s;
	// curandStateMRG32k3a s;
	curand_init(sbi, 0, 0, &s);

	float x = curand_uniform(&s);
	tau[sbi] = x;
}

void h_fill_tau_array(thrust::device_vector<float> &tau)
{
	fill_tau_array<<<1, tau.size()>>>(thrust::raw_pointer_cast(tau.data()), tau.size());
}

int h_get_min_tau(thrust::device_vector<float> &tau)
{
	thrust::device_vector<float>::iterator iter = thrust::min_element(tau.begin(), tau.end());
	return iter - tau.begin();
}
