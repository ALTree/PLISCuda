#include "../include/cuda/ssa.cuh"

__device__ int choose_rand_reaction(rates rates, float rand)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return -1;

	// if R (the sum of the reaction rates) is zero,
	// we can't fire any reaction
	if (rates.matrix[GET_RATE(0, sbi)] == 0)
		return -1;

	float sum = rates.matrix[sbi];
	float scaled_sum = sum * rand;
	float partial_sum = 0;

	int ri = 0;
	while (partial_sum <= scaled_sum) {
		partial_sum += rates.reaction[GET_RR(ri, sbi)];
		ri++;
	}

	return ri - 1;
}

__device__ int choose_rand_specie(neigh neigh, rates rates, float rand)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return -1;

	// if D (the sum of the diffusion rates) is zero,
	// we can't diffuse any specie
	if (rates.matrix[GET_RATE(1, sbi)] == 0)
		return -1;

	int neigh_count = neigh.count[sbi];

	// we need to scale back rate_matrix[2][sbi] before performing
	// the linear scaling
	float sum = rates.matrix[SBC * 1 + sbi] / neigh_count;
	float scaled_sum = sum * rand;
	float partial_sum = 0;

	int spi = 0;
	while (partial_sum <= scaled_sum) {
		partial_sum += rates.diffusion[GET_DR(spi, sbi)];
		spi++;
	}

	return spi - 1;
}

__global__ void initialize_prngstate_array(curandStateMRG32k3a * prngstate)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

#ifndef TEST
	curand_init(clock64() * sbi, 0, 0, &prngstate[sbi]);
#else 
	curand_init(sbi, 0, 0, &prngstate[sbi]);
#endif

}

__global__ void ssa_step(state state, reactions reactions, neigh neigh, rates rates,
						 int min_sbi, float * current_time, char * leap, curandStateMRG32k3a * s)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC || leap[sbi] == LEAP_CR || leap[sbi] == LEAP_NOCR || min_sbi != sbi)
		return;

	float rand = curand_uniform(&s[sbi]);

	if (rand < rates.matrix[GET_RATE(0, sbi)] / rates.matrix[GET_RATE(2, sbi)]) {    // reaction

		int ri = choose_rand_reaction(rates, rand);
#ifdef LOG
		printf("(%f) [subv %d] fire reaction %d  [SSA]\n", *current_time, sbi, ri);
#endif
			
		// fire reaction and update the state of the system
		// if sbi = min_sbi then it should be guaranteed that ri != -1
		for (int spi = 0; spi < SPC; spi++) {
			int delta = reactions.p[GET_COEFF(spi, ri)] - reactions.r[GET_COEFF(spi, ri)];
			atomicAdd(&state.next[GET_SPI(spi, sbi)], delta);
		}

		// set our own OP to SSA (we can't fast-forward if we fired a reaction)
		leap[sbi] = SSA;

	} else {    // diffusion

		int spi = choose_rand_specie(neigh, rates, rand);

		// choose a random destination
		int rdi;
		do {
			rdi = (int) (curand_uniform(&s[sbi]) * 6);
		} while (rdi > 5);

		// get index of neighbour #rdi (overwrite rdi, whatever)
		rdi = neigh.index[sbi * 6 + rdi];

#ifdef LOG
		printf("(%f) [subv %d] diffuse specie %d in subvolume %d  [SSA]\n", *current_time, sbi, spi, rdi);
#endif
			
		// If rdi == sbi (i.e. diffuse to myself) don't do anything
		if (rdi != sbi) {
			atomicSub(&state.next[GET_SPI(spi, sbi)], 1);
			atomicAdd(&state.next[GET_SPI(spi, rdi)], 1);
			if (leap[rdi] == SSA_FF)
				leap[rdi] = SSA;    // set the OP of the receiver to SSA
		}

	}

}

