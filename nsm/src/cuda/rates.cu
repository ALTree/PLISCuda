#include "../../include/cuda/rates.cuh"

// #define DEBUG

// TODO: on multiple threads
__device__ float react_rate(int * state, int * reactants, int sbc, int spc, int rc, int sbi, int ri, float * rrc)
{
#ifdef DEBUG
	printf("---------- begin react_rate( ) ---------- \n");
	printf("#reactions = %d, #species = %d, #subs = %d\n", rc, spc, sbc);
	printf("sub_index = %d, reaction_index = %d\n", sbi, ri);
#endif

	// search for the first specie in the reactions array that
	// does have a positive coefficent
	int index1 = ri;
	int specie_index = 0;
	while (reactants[index1] == 0) {
		index1 += rc;
		specie_index++;
	}

	if (reactants[index1] == 2) {    // bi_same reaction type
		// get specie count for that specie in the current subvolume
		// int specie_count = state[specie_index * subvolumes_count + subvolume_index];
		int specie_count = state[CUDA_GET_SPI(specie_index, sbi, sbc)];
#ifdef DEBUG
		printf("hit bi_same. Specie_index = %d, specie_count = %d\n", specie_index, species_count);
		printf("----------   end react_rate( ) ---------- \n\n");
#endif
		return 0.5 * specie_count * (specie_count - 1) * rrc[ri];
	}

	// if specie_index == # of species we are in a uni reaction
	if (specie_index != spc - 1) {

		// search for a possibile other specie with positive coefficient
		int index2 = index1 + rc;
		int specie_index2 = specie_index + 1;
		while (reactants[index2] == 0 && index2 < spc * rc) {
			index2 += rc;
			specie_index2++;
		}

		if (reactants[index2] != 0) {    // bi_diff reaction type
			int specie1_count = state[CUDA_GET_SPI(specie_index, sbi, sbc)];
			int specie2_count = state[CUDA_GET_SPI(specie_index2, sbi, sbc)];
#ifdef DEBUG
			printf("hit bi_diff. Specie_index1 = %d, specie_index2 = %d,\n", specie_index, specie_index2);
			printf("    specie1_count = %d, specie2_count = %d\n", specie1_count, specie2_count);
			printf("----------   end react_rate( ) ----------\n\n");
#endif
			return specie1_count * specie2_count * rrc[ri];
		}
	}

	// uni reaction type
	int specie_count = state[CUDA_GET_SPI(specie_index, sbi, sbc)];
#ifdef DEBUG
	printf("hit uni. Specie_index = %d, specie_count = %d, ", specie_index, specie_count);
	printf("----------   end react_rate( ) ---------- \n\n");
#endif
	return specie_count * rrc[ri];
}

__device__ void react_rates(int * state, int * reactants, int sbc, int spc, int rc, float * rrc,
		float * react_rates_array)
{
	int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbc <= sbi) {
		return;
	}

	for (int i = 0; i < rc; i++) {
		react_rates_array[sbc * i + sbi] = react_rate(state, reactants, sbc, spc, rc, sbi, i, rrc);
	}

}

// TODO: on multiple threads
__device__ void diff_rates(int * state, int sbc, int spc, int sbi, float * drc, float * result)
{
	for (int i = 0; i < spc; i++) {
		result[i] = drc[i] * state[CUDA_GET_SPI(i, sbi, sbc)];
	}

}
