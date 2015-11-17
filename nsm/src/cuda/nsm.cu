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

__device__ int choose_rand_specie(int * topology, float * rate_matrix, float * diff_rates_array, float rand)
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

__global__ void fill_tau_array(float * tau, float * rate_matrix)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	curandStateMRG32k3a s;
	curand_init(2 * sbi, 0, 0, &s);    // initialize with *2sbi to avoid getting the same first value
									   // later when we use curand_init(sbi, ..)

	float rand = curand_uniform(&s);
	tau[sbi] = -logf(rand) / rate_matrix[GET_RATE(2, sbi)];
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

__global__ void nsm_step(int * state, int * reactants, int * products, int * topology, float * rate_matrix, float * rrc,
		float * drc, float * react_rates_array, float * diff_rates_array, float * tau, int min_sbi,
		curandStateMRG32k3a * s)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	float rand = curand_uniform(&s[sbi]);

#ifdef DEBUG
	printf("[sbv %d] tau = %f, rand = %f\n", sbi, tau[sbi], rand);
#endif

	if (rand < rate_matrix[GET_RATE(0, sbi)] / rate_matrix[GET_RATE(2, sbi)]) {    // reaction

	// choose a random reaction to fire
		int ri = choose_rand_reaction(rate_matrix, react_rates_array, rand);

#ifdef DEBUG
		if (ri >= RC) {
			printf(">>>>> ARGH! @ [subv %d]: random reaction index = %d\n", sbi, ri);
		}
#endif

		// ri = -1 means we can't fire any reaction in this subvolume
		if (sbi == min_sbi && ri != -1) {
			printf("(%f) [subv %d] fire reaction %d\n", tau[sbi], sbi, ri);
		}

		// fire reaction and update the state of the system
		// if sbi = min_sbi then it should be guaranteed that ri != -1
		// TODO: check(?)
		if (sbi == min_sbi) {    // (but only if you are the choosen one)
			for (int i = 0; i < SPC; i++)
				state[GET_SPI(i, sbi)] += products[i * RC + ri] - reactants[i * RC + ri];
		}

		// TODO: do we need this?
		__syncthreads();

		// update rate matrix
		react_rates(state, reactants, rrc, react_rates_array);
		diff_rates(state, drc, diff_rates_array);
		update_rate_matrix(topology, rate_matrix, react_rates_array, diff_rates_array);
	} else {    // diffusion

		// choose a random specie to diffuse
		int spi = choose_rand_specie(topology, rate_matrix, diff_rates_array, rand);

#ifdef DEBUG
		if (spi >= SPC) {
			printf(">>>>> ARGH! @ [subv %d] random specie index = %d\n", sbi, spi);
		}
#endif

		// choose a random destination
		// TODO: we need to re-use the rand we already have.
		int rdi;
		do {
			rdi = (int) (curand_uniform(&s[sbi]) * 6);
		} while (rdi > 5);

		// get index of neighbour #rdi (overwrite rdi, whatever)
		rdi = topology[sbi * 6 + rdi];

#ifdef DEBUG
		if (rdi >= SBC) {
			printf(">>>>> ARGH! @ [subv %d] random neigh = %d\n", sbi, rdi);
		}
#endif

		if (sbi == min_sbi) {
			printf("(%f) [subv %d] diffuse specie %d in subvolume %d\n", tau[sbi], sbi, spi, rdi);
		}

		// Update state iff we are the choosen one.
		// Also if we rdi == sbi (i.e. diffuse to myself) don't do anything
		if (sbi == min_sbi && rdi != sbi) {
			state[GET_SPI(spi, sbi)] -= 1;
			state[GET_SPI(spi, rdi)] += 1;
		}

		// TODO: do we need this?
		__syncthreads();

		// update rate matrix
		// The destination subvolume will update its own rates... right?
		react_rates(state, reactants, rrc, react_rates_array);
		diff_rates(state, drc, diff_rates_array);
		update_rate_matrix(topology, rate_matrix, react_rates_array, diff_rates_array);
	}

	// compute next event time for this subvolume
	rand = curand_uniform(&s[sbi]);
	tau[sbi] = -logf(rand) / rate_matrix[GET_RATE(2, sbi)] + tau[min_sbi];
}
