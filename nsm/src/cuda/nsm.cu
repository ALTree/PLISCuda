#include "../../include/cuda/nsm.cuh"

// #define DEBUG

__device__ int choose_rand_reaction(float * rate_matrix, float * react_rates_array, float rand)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return -1;

	// if R (the sum of the reaction rates) is zero,
	// we can't fire any reaction
	if (rate_matrix[GET_RATE(0, sbi)] == 0)
		return -1;

	float sum = rate_matrix[sbi];
	float scaled_sum = sum * rand;
	float partial_sum = 0;

	int ri = 0;
	while (partial_sum <= scaled_sum) {
		partial_sum += react_rates_array[ri * SBC + sbi];
		ri++;
	}

	return ri - 1;
}

__device__ int choose_rand_specie(unsigned int * topology, float * rate_matrix, float * diff_rates_array, float rand)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return -1;

	// if D (the sum of the diffusion rates) is zero,
	// we can't diffuse any specie
	if (rate_matrix[GET_RATE(1, sbi)] == 0)
		return -1;

	int neigh_count = 0;
	for (int i = 0; i < 6; i++)
		neigh_count += (topology[sbi * 6 + i] != sbi);

	// we need to scale back rate_matrix[2][sbi] before performing
	// the linear scaling

	float sum = rate_matrix[SBC * 1 + sbi] / neigh_count;
	float scaled_sum = sum * rand;
	float partial_sum = 0;

	int spi = 0;
	while (partial_sum <= scaled_sum) {
		partial_sum += diff_rates_array[spi * SBC + sbi];
		spi++;
	}

	return spi - 1;
}

// TODO: remove fill_tau_array and use this for everything
__global__ void fill_prngstate_array(curandStateMRG32k3a * prngstate)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	curand_init(sbi, 0, 0, &prngstate[sbi]);
}

int h_get_min_tau(thrust::device_vector<float> &tau)
{
	thrust::device_vector<float>::iterator iter = thrust::min_element(tau.begin(), tau.end());
	return iter - tau.begin();
}

__global__ void nsm_step(int * state, int * reactants, int * products, unsigned int * topology, float * rate_matrix,
		float * rrc, float * drc, float * react_rates_array, float * diff_rates_array, float * tau, int min_sbi,
		float * current_time, bool * leap, curandStateMRG32k3a * s)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC || leap[sbi] || min_sbi != sbi)
		return;

	float rand = curand_uniform(&s[sbi]);

	if (rand < rate_matrix[GET_RATE(0, sbi)] / rate_matrix[GET_RATE(2, sbi)]) {    // reaction

	// choose a random reaction to fire
		int ri = choose_rand_reaction(rate_matrix, react_rates_array, rand);

		// ri = -1 means we can't fire any reaction in this subvolume
		if (ri != -1) {
			printf("(%f) [subv %d][SSA] fire reaction %d\n", *current_time, sbi, ri);
		} else {    // TODO: remove
			printf("sbv %d], got ri == -1, rand = %f, rate[0, sbi] = %f\n", sbi, rand, rate_matrix[GET_RATE(0, sbi)]);
		}

		// fire reaction and update the state of the system
		// if sbi = min_sbi then it should be guaranteed that ri != -1
		// TODO: check(?)
		if (ri != -1) {
			for (int i = 0; i < SPC; i++)
				state[GET_SPI(i, sbi)] += products[GET_COEFF(i, ri)] - reactants[GET_COEFF(i, ri)];
		}

	} else {    // diffusion

		// choose a random specie to diffuse
		int spi = choose_rand_specie(topology, rate_matrix, diff_rates_array, rand);

		// choose a random destination
		// TODO: we need to re-use the rand we already have.
		int rdi;
		do {
			rdi = (int) (curand_uniform(&s[sbi]) * 6);
		} while (rdi > 5);

		// get index of neighbour #rdi (overwrite rdi, whatever)
		rdi = topology[sbi * 6 + rdi];

		printf("(%f) [subv %d][SSA] diffuse specie %d in subvolume %d\n", *current_time, sbi, spi, rdi);

		// Update state iff we are the choosen one.
		// Also if we rdi == sbi (i.e. diffuse to myself) don't do anything
		if (rdi != sbi) {
			state[GET_SPI(spi, sbi)] -= 1;
			state[GET_SPI(spi, rdi)] += 1;
		}

	}

}
