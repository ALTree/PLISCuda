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
	while (partial_sum <= scaled_sum) {
		partial_sum += react_rates_array[ri * sbc + sbi];
		ri++;
	}

	return ri - 1;
}

__device__ int choose_rand_specie(int sbc, int spc, float * rate_matrix, float * diff_rates_array, float rand)
{
	int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= sbc)
		return -1;

	float sum = rate_matrix[sbc * 1 + sbi];
	float scaled_sum = sum * rand;
	float partial_sum = 0;

	int spi = 0;
	while (partial_sum <= scaled_sum) {
		partial_sum += diff_rates_array[spi * sbc + sbi];
		spi++;
	}

	return spi - 1;
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

int h_get_min_tau(thrust::device_vector<float> &tau)
{
	thrust::device_vector<float>::iterator iter = thrust::min_element(tau.begin(), tau.end());
	return iter - tau.begin();
}

__global__ void nsm_step(int * state, int * reactants, int * products, int * topology, int sbc, int spc, int rc,
		float * rate_matrix, float * react_rates_array, float * diff_rates_array, float * tau, int min_sbi)
{
	int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= sbc)
		return;

	// create and initialize thread's prng
	curandStateMRG32k3a s;
	curand_init(sbi, 0, 0, &s);

	float rand = curand_uniform(&s);

	if (rand < rate_matrix[GET_RATE(0, sbi, sbc)] / rate_matrix[GET_RATE(2, sbi, sbc)]) {
		// fire a reaction
		printf("(%f) [subvolume %d] fire random reaction\n", tau[sbi], sbi);

		// choose a random reaction to fire
		int ri = choose_rand_reaction(sbc, rc, rate_matrix, react_rates_array, rand);
		if (ri < 0 || ri >= rc) {
			printf(">>>>>>>>>>>>>>>> ARGH! @ sbi: random reaction index = %d\n", sbi, ri);
		}

		// fire reaction and update the state of the system
		if (sbi == min_sbi) {    // (but only if you are the choosen one)
			for (int i = 0; i < spc; i++)
				state[GET_SPI(i, sbi, sbc)] += products[i * rc + ri] - reactants[i * rc + ri];
		}

		// update rate matrix
		update_rate_matrix(topology, sbc, spc, rc, rate_matrix, react_rates_array, diff_rates_array);

		// compute next event time for this subvolume
		rand = curand_uniform(&s);
		tau[sbi] += -logf(rand)/rate_matrix[GET_RATE(2, sbi, sbc)];


	} else {
		// diffuse a specie
		printf("(%f) [subvolume %d] diffuse random specie (not implemented)\n", tau[sbi], sbi);
		rand = curand_uniform(&s);
		tau[sbi] += -logf(rand)/rate_matrix[GET_RATE(2, sbi, sbc)];
	}
}
