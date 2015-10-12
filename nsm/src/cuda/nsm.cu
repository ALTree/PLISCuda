#include "../../include/cuda/nsm.cuh"

// #define DEBUG

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

__device__ float * react_rates(int * state, int * reactants, int sbc, int spc, int rc, int sbi, float * rrc)
{
	__shared__ extern float react_rates_array[];    // we need extern because the size is not a compile-time constant
													// we'll need to allocate during the kernel invocation
													// TODO: rethink about this

	for (int i = 0; i < rc; i++) {
		react_rates_array[i] = react_rate(state, reactants, sbc, spc, rc, sbi, i, rrc);
	}

	return react_rates_array;
}

__device__ float * diff_rates(int * state, int subvolumes_count, int species_count, int subvolume_index,
		float * diffusion_rates_constants)
{
	__shared__ extern float diffusion_rates_array[];
	for (int i = 0; i < species_count; i++) {
		diffusion_rates_array[i] = diffusion_rates_constants[i]
				* state[CUDA_GET_SPI(i, subvolume_index, subvolumes_count)];
	}

	return diffusion_rates_array;
}

/*
__device__ void rate_matrix_row(int * state, int * reactants, int subvolumes_count, int species_count,
		int reactions_count, float * reaction_rate_constants, float * diffusion_rate_constants, float * rate_matrix,
		int subvolume_index)
{
	// compute new reaction rates
	float * react_rates_array = react_rates(reactants, reactions_count, state, subvolumes_count, species_count,
			subvolume_index, reaction_rate_constants);
	float reactions_rates_sum = sum_fp_array(react_rates_array, reactions_count);

	// compute new diffusion rates
	float * diff_rates_array = diff_rates(state, subvolumes_count, species_count, subvolume_index,
			diffusion_rate_constants);
	float diffusion_rates_sum = sum_fp_array(diff_rates_array, species_count);

	// update rate matrix
	rate_matrix[subvolume_index] = reactions_rates_sum;
	rate_matrix[subvolume_index * 2] = diffusion_rates_sum;
	rate_matrix[subvolume_index * 3] = reactions_rates_sum + diffusion_rates_sum;
}
*/
