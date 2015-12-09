#include "../../include/cuda/ssa.cuh"

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
		partial_sum += react_rates_array[GET_RR(ri, sbi)];
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
		partial_sum += diff_rates_array[GET_DR(spi, sbi)];
		spi++;
	}

	return spi - 1;
}

// TODO: remove fill_tau_array and use this for everything
__global__ void initialize_prngstate_array(curandStateMRG32k3a * prngstate)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	curand_init(sbi, 0, 0, &prngstate[sbi]);
}

__global__ void ssa_step(int * state, int * reactants, int * products, unsigned int * topology, float * rate_matrix,
		float * react_rates_array, float * diff_rates_array, int min_sbi, float * current_time, bool * leap,
		curandStateMRG32k3a * s)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC || leap[sbi] || min_sbi != sbi)
		return;

	float rand = curand_uniform(&s[sbi]);

	if (rand < rate_matrix[GET_RATE(0, sbi)] / rate_matrix[GET_RATE(2, sbi)]) {    // reaction

		int ri = choose_rand_reaction(rate_matrix, react_rates_array, rand);
		printf("(%f) [subv %d][SSA] fire reaction %d\n", *current_time, sbi, ri);

		// fire reaction and update the state of the system
		// if sbi = min_sbi then it should be guaranteed that ri != -1
		for (int i = 0; i < SPC; i++)
			state[GET_SPI(i, sbi)] += products[GET_COEFF(i, ri)] - reactants[GET_COEFF(i, ri)];

	} else {    // diffusion

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

		// If rdi == sbi (i.e. diffuse to myself) don't do anything
		// TODO: atomic?
		if (rdi != sbi) {
			state[GET_SPI(spi, sbi)] -= 1;
			state[GET_SPI(spi, rdi)] += 1;
		}

	}

}

