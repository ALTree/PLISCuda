#include "../../include/cuda/nsm.cuh"

#define DEBUG

__device__ int choose_rand_reaction(int sbc, int rc, float * rate_matrix, float * react_rates_array, float rand)
{
	int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= sbc)
		return -1;

	// if R (the sum of the reaction rates) is zero,
	// we can't fire any reaction
	if (rate_matrix[GET_RATE(0, sbi, sbc)] == 0)
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

__device__ int choose_rand_specie(int * topology, int sbc, int spc, float * rate_matrix, float * diff_rates_array,
		float rand)
{
	int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= sbc)
		return -1;

	// if D (the sum of the diffusion rates) is zero,
	// we can't diffuse any specie
	if (rate_matrix[GET_RATE(1, sbi, sbc)] == 0)
		return -1;

	int neigh_count = 0;
	for (int i = 0; i < 6; i++)
		neigh_count += (topology[sbi * 6 + i] != -1);

	// we need to scale back rate_matrix[2][sbi] before performing
	// the linear scaling

	float sum = rate_matrix[sbc * 1 + sbi] / neigh_count;
	float scaled_sum = sum * rand;
	float partial_sum = 0;

	int spi = 0;
	while (partial_sum <= scaled_sum) {
		partial_sum += diff_rates_array[spi * sbc + sbi];
		spi++;
	}

	return spi - 1;
}

__global__ void fill_tau_array(float * tau, float * rate_matrix, int sbc)
{
	int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= sbc)
		return;

	curandStateMRG32k3a s;
	curand_init(2 * sbi, 0, 0, &s);    // initialize with *2sbi to avoid getting the same first value
									   // later when we use curand_init(sbi, ..)

	float rand = curand_uniform(&s);
	tau[sbi] = -logf(rand) / rate_matrix[GET_RATE(2, sbi, sbc)];
}

int h_get_min_tau(thrust::device_vector<float> &tau)
{
	thrust::device_vector<float>::iterator iter = thrust::min_element(tau.begin(), tau.end());
	return iter - tau.begin();
}

__global__ void nsm_step(int * state, int * reactants, int * products, int * topology, int sbc, int spc, int rc,
		float * rate_matrix, float * rrc, float * drc, float * react_rates_array, float * diff_rates_array, float * tau,
		int min_sbi, int step)
{
	int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= sbc)
		return;

	// create and initialize thread's prng
	curandStateMRG32k3a s;
	curand_init(sbi, 0, step, &s);

	float rand = curand_uniform(&s);

#ifdef DEBUG
	printf("[sbv %d] tau = %f, rand = %f\n", sbi, tau[sbi], rand);
#endif

	if (rand < rate_matrix[GET_RATE(0, sbi, sbc)] / rate_matrix[GET_RATE(2, sbi, sbc)]) {
		// fire a reaction

		// choose a random reaction to fire
		int ri = choose_rand_reaction(sbc, rc, rate_matrix, react_rates_array, rand);

		if (ri == -1)    // we can't fire any reaction in this subvolume
			goto UPDATE_TAU;

		if (ri >= rc) {
			printf(">>>>>>>>>>>>>>>> ARGH! @ [subv %d]: random reaction index = %d\n", sbi, ri);
		}

		if (sbi == min_sbi) {
			printf("(%f) [subv %d] fire reaction %d\n", tau[sbi], sbi, ri);
		}

		// fire reaction and update the state of the system
		if (sbi == min_sbi) {    // (but only if you are the choosen one)
			for (int i = 0; i < spc; i++)
				state[GET_SPI(i, sbi, sbc)] += products[i * rc + ri] - reactants[i * rc + ri];
		}

		// TODO: do we need this?
		__syncthreads();

		// update rate matrix
		react_rates(state, reactants, sbc, spc, rc, rrc, react_rates_array);
		diff_rates(state, sbc, spc, drc, diff_rates_array);
		update_rate_matrix(topology, sbc, spc, rc, rate_matrix, react_rates_array, diff_rates_array);
	} else {
		// diffuse a specie

		// choose a random specie to diffuse
		int spi = choose_rand_specie(topology, sbc, spc, rate_matrix, diff_rates_array, rand);

		if (spi == -1)    // we can't diffuse any specie in this subvolume
			goto UPDATE_TAU;

		if (spi >= spc) {
			printf(">>>>>>>>>>>>>>>> ARGH! @ [subv %d: random specie index = %d\n", sbi, spi);
		}

		// choose a random destination
		// TODO: we need to re-use the rand we already have.
		//       Also find a better way to ensure fairness on
		//       index 5.
		int rdi = (int) curand_uniform(&s) * 6;
		while (rdi > 5)
			rdi = (int) curand_uniform(&s);

		// get index of neighbour #rdi (overwrite rdi, whatever)
		rdi = topology[sbi * 6 + rdi];

		if (rdi < 0 || rdi >= sbc) {
			printf(">>>>>>>>>>>>>>>> ARGH! @ sbi: random neigh index = %d\n", sbi, rdi);
		}

		if (sbi == min_sbi) {
			printf("(%f) [subv %d] diffuse specie %d in subvolume %d\n", tau[sbi], sbi, spi, rdi);
		}

		// Update state iff we are the choosen one.
		// Also if we hit a -1 (i.e. diffuse to myself) don't do anything
		if (sbi == min_sbi && rdi != -1) {
			state[GET_SPI(spi, sbi, sbc)] -= 1;
			state[GET_SPI(spi, rdi, sbc)] += 1;
		}

		// TODO: do we need this?
		__syncthreads();

		// update rate matrix
		// The destination subvolume will update its own rates... right?
		react_rates(state, reactants, sbc, spc, rc, rrc, react_rates_array);
		diff_rates(state, sbc, spc, drc, diff_rates_array);
		update_rate_matrix(topology, sbc, spc, rc, rate_matrix, react_rates_array, diff_rates_array);
	}

	// compute next event time for this subvolume
	UPDATE_TAU: rand = curand_uniform(&s);
	tau[sbi] = -logf(rand) / rate_matrix[GET_RATE(2, sbi, sbc)] + tau[min_sbi];
}
