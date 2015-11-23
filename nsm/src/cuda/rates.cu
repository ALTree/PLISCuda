#include "../../include/cuda/rates.cuh"

__device__ float react_rate(int * state, int * reactants, int sbi, int ri, float * rrc)
{
	// search for the first specie in the reactions array that
	// does have a positive coefficent
	int index1 = ri;
	int specie_index = 0;
	while (reactants[index1] == 0) {
		index1 += RC;
		specie_index++;
	}

	if (reactants[index1] == 2) {    // bi_same reaction type
		// get specie count for that specie in the current subvolume
		// int specie_count = state[specie_index * subvolumes_count + subvolume_index];
		int specie_count = state[GET_SPI(specie_index, sbi)];
		return 0.5 * specie_count * (specie_count - 1) * rrc[ri];
	}

	// if specie_index == # of species we are in a uni reaction
	if (specie_index != SPC - 1) {

		// search for a possibile other specie with positive coefficient
		int index2 = index1 + RC;
		int specie_index2 = specie_index + 1;
		while (reactants[index2] == 0 && index2 < (SPC * RC - 1)) {
			index2 += RC;
			specie_index2++;
		}

		if (reactants[index2] != 0) {    // bi_diff reaction type
			int specie1_count = state[GET_SPI(specie_index, sbi)];
			int specie2_count = state[GET_SPI(specie_index2, sbi)];
			return specie1_count * specie2_count * rrc[ri];
		}
	}

	// uni reaction type
	int specie_count = state[GET_SPI(specie_index, sbi)];
	return specie_count * rrc[ri];
}

__device__ void react_rates(int * state, int * reactants, float * rrc, float * react_rates_array)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	for (int i = 0; i < RC; i++) {
		react_rates_array[SBC * i + sbi] = react_rate(state, reactants, sbi, i, rrc);
	}
}

__device__ void diff_rates(int * state, float * drc, float * diff_rates_array)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	for (int i = 0; i < SPC; i++) {
		diff_rates_array[SBC * i + sbi] = drc[i] * state[GET_SPI(i, sbi)];
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
	for (int i = 0; i < RC; i++)
		react_sum += react_rates_array[SBC * i + sbi];

	// sum diffusion rates
	float diff_sum = 0.0;
	for (int i = 0; i < SPC; i++)
		diff_sum += diff_rates_array[SBC * i + sbi];

	// count subvolume neighbours (since diff_rate = #neighbours x diff_sum)
	// TODO: write # of neighbours somewhere and use it so we can remove the
	//       topology parameter
	int neigh_count = 0;
	for (int i = 0; i < 6; i++)
		neigh_count += (topology[sbi * 6 + i] != sbi);

	diff_sum *= neigh_count;

	// write data into rate matrix
	rate_matrix[SBC * 0 + sbi] = react_sum;
	rate_matrix[SBC * 1 + sbi] = diff_sum;
	rate_matrix[SBC * 2 + sbi] = react_sum + diff_sum;
}

__global__ void compute_rates(int * state, int * reactants, unsigned int * topology, float * rate_matrix, float * rrc,
		float * drc, float * react_rates_array, float * diff_rates_array)
{
	react_rates(state, reactants, rrc, react_rates_array);
	diff_rates(state, drc, diff_rates_array);
	update_rate_matrix(topology, rate_matrix, react_rates_array, diff_rates_array);
}
