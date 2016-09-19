#include "../include/cuda/rates.cuh"

__device__ float react_rate(state state, reactions reactions, int sbi, int ri, rates rates)
{
	// search for the first specie in the reactions array that
	// does have a positive coefficent
	int i = ri;
	int spi1 = 0;
	while (reactions.r[i] == 0 && spi1 < SPC) {
		i += RC;
		spi1++;
	}

	if (reactions.r[i] == 2) {    // bi_same reaction type
		int sp_count = state.curr[GET_SPI(spi1, sbi)];
		return 0.5 * (rates.rc[ri] *  sp_count) * (sp_count - 1); // careful with overflow
	}

	// if we didn't look at all the species yet, search
	// for a possible second positive value
	if (spi1 != SPC - 1) {
		int j = i + RC; // start from the next specie
		int spi2 = spi1 + 1;
		while (spi2+1 < SPC && reactions.r[j] == 0) {
			j += RC;
			spi2++;
		}

		if (reactions.r[j] != 0) {    // bi_diff reaction type
			int sp1_count = state.curr[GET_SPI(spi1, sbi)];
			int sp2_count = state.curr[GET_SPI(spi2, sbi)];
			return (rates.rc[ri] * sp1_count) * sp2_count; // careful with overflow
		}
	}
	// uni reaction type

	int sp_count = state.curr[GET_SPI(spi1, sbi)];
	return sp_count * rates.rc[ri];
}

__device__ void react_rates(state state, reactions reactions, rates rates)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	for (int ri = 0; ri < RC; ri++) {
		rates.reaction[GET_RR(ri, sbi)] = react_rate(state, reactions, sbi, ri, rates);
	}
}

__device__ void diff_rates(state state, rates rates)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	for (int spi = 0; spi < SPC; spi++) {
		rates.diffusion[GET_DR(spi, sbi)] = rates.dc[spi] * state.curr[GET_SPI(spi, sbi)];
	}
}

__device__ void update_rate_matrix(unsigned int * topology, rates rates)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	// sum reaction rates
	float react_sum = 0.0;
	for (int ri = 0; ri < RC; ri++)
		react_sum += rates.reaction[GET_RR(ri, sbi)];

	// sum diffusion rates
	float diff_sum = 0.0;
	for (int spi = 0; spi < SPC; spi++)
		diff_sum += rates.diffusion[GET_DR(spi, sbi)];

	// count subvolume neighbours (since diff_rate = #neighbours x diff_sum)
	// TODO: write # of neighbours somewhere and use it so we can remove the
	//       topology parameter
	int neigh_count = 0;
	for (int i = 0; i < 6; i++)
		neigh_count += (topology[sbi * 6 + i] != sbi);

	diff_sum *= neigh_count;

	rates.matrix[GET_RATE(0, sbi)] = react_sum;
	rates.matrix[GET_RATE(1, sbi)] = diff_sum;
	rates.matrix[GET_RATE(2, sbi)] = react_sum + diff_sum;
}

__global__ void compute_rates(state state, reactions reactions, unsigned int * topology, rates rates,
							  int * d_subv_consts)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	// find the position of the constants set we
	// have to use (bc compartimentation support)
	int dsc = d_subv_consts[sbi];
	rates.rc = &rates.rc[dsc*RC];
	rates.dc = &rates.dc[dsc*SPC];

	react_rates(state, reactions, rates);
	diff_rates(state, rates);
	update_rate_matrix(topology, rates);
}
