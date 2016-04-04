#include "../include/cuda/rates.cuh"

__device__ float react_rate(int * state, int * reactants, int sbi, int ri, float * rrc)
{
	// search for the first specie in the reactions array that
	// does have a positive coefficent
	int i = ri;
	int spi1 = 0;
	while (reactants[i] == 0) {
		i += RC;
		spi1++;
	}

	if (reactants[i] == 2) {    // bi_same reaction type
		int sp_count = state[GET_SPI(spi1, sbi)];
		return 0.5 * sp_count * (sp_count - 1) * rrc[ri];
	}

	// if we didn't look at all the species yet, search
	// for a possible second positive value
	if (spi1 != SPC - 1) {
		int j = i + RC; // start from the next specie
		int spi2 = spi1 + 1;
		while (reactants[j] == 0 && j < (SPC * RC - 1)) {
			j += RC;
			spi2++;
		}

		if (reactants[j] != 0) {    // bi_diff reaction type
			int sp1_count = state[GET_SPI(spi1, sbi)];
			int sp2_count = state[GET_SPI(spi2, sbi)];
			return (rrc[ri] * sp1_count) * sp2_count; // careful with overflow
		}
	}

	// uni reaction type
	int sp_count = state[GET_SPI(spi1, sbi)];
	return sp_count * rrc[ri];
}

__device__ void react_rates(int * state, int * reactants, float * rrc, float * react_rates_array)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	for (int ri = 0; ri < RC; ri++) {
		react_rates_array[GET_RR(ri, sbi)] = react_rate(state, reactants, sbi, ri, rrc);
	}
}

__device__ void diff_rates(int * state, float * drc, float * diff_rates_array)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	for (int spi = 0; spi < SPC; spi++) {
		diff_rates_array[GET_DR(spi, sbi)] = drc[spi] * state[GET_SPI(spi, sbi)];
	}
}

__device__ void update_rate_matrix(unsigned int * topology, float * rate_matrix, float * react_rates_array,
		float * diff_rates_array)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	// sum reaction rates
	float react_sum = 0.0;
	for (int ri = 0; ri < RC; ri++)
		react_sum += react_rates_array[GET_RR(ri, sbi)];

	// sum diffusion rates
	float diff_sum = 0.0;
	for (int spi = 0; spi < SPC; spi++)
		diff_sum += diff_rates_array[GET_DR(spi, sbi)];

	// count subvolume neighbours (since diff_rate = #neighbours x diff_sum)
	// TODO: write # of neighbours somewhere and use it so we can remove the
	//       topology parameter
	int neigh_count = 0;
	for (int i = 0; i < 6; i++)
		neigh_count += (topology[sbi * 6 + i] != sbi);

	diff_sum *= neigh_count;

	rate_matrix[GET_RATE(0, sbi)] = react_sum;
	rate_matrix[GET_RATE(1, sbi)] = diff_sum;
	rate_matrix[GET_RATE(2, sbi)] = react_sum + diff_sum;
}

__global__ void compute_rates(int * state, int * reactants, unsigned int * topology, float * rate_matrix, float * rrc,
		float * drc, int * d_subv_consts, float * react_rates_array, float * diff_rates_array)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	float * rrcp = &rrc[d_subv_consts[sbi]*RC];
	float * drcp = &drc[d_subv_consts[sbi]*SPC];

	react_rates(state, reactants, rrcp, react_rates_array);
	diff_rates(state, drcp, diff_rates_array);
	update_rate_matrix(topology, rate_matrix, react_rates_array, diff_rates_array);
}
