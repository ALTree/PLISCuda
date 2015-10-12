#include "../../include/cuda/nsm.cuh"

// #define DEBUG

__device__ float react_rate(int * reactants, int reactions_count, int * state, int species_count, int subvolumes_count,
		int subvolume_index, float * reaction_rate_constants, int reaction_number)
{
#ifdef DEBUG
	printf("---------- begin react_rate( ) ---------- \n");
	printf("#reactions = %d, #species = %d, #subs = %d\n", reactions_count, species_count, subvolumes_count);
	printf("sub_index = %d, reaction_number = %d\n", subvolume_index, reaction_number);
#endif

	// search for the first specie in the reactions array that
	// does have a positive coefficent
	int index1 = reaction_number;
	int specie_index = 0;
	while (reactants[index1] == 0) {
		index1 += reactions_count;
		specie_index++;
	}

	if (reactants[index1] == 2) {    // bi_same reaction type
		// get specie count for that specie in the current subvolume
		// int specie_count = state[specie_index * subvolumes_count + subvolume_index];
		int specie_count = state[CUDA_GET_SPI(specie_index, subvolume_index, subvolumes_count)];
#ifdef DEBUG
		printf("hit bi_same. Specie_index = %d, specie_count = %d\n", specie_index, species_count);
		printf("----------   end react_rate( ) ---------- \n\n");
#endif
		return 0.5 * specie_count * (specie_count - 1) * reaction_rate_constants[reaction_number];
	}

	// if specie_index == # of species we are in a uni reaction
	if (specie_index != species_count - 1) {

		// search for a possibile other specie with positive coefficient
		int index2 = index1 + reactions_count;
		int specie_index2 = specie_index + 1;
		while (reactants[index2] == 0 && index2 < species_count * reactions_count) {
			index2 += reactions_count;
			specie_index2++;
		}

		if (reactants[index2] != 0) {    // bi_diff reaction type
			int specie1_count = state[CUDA_GET_SPI(specie_index, subvolume_index, subvolumes_count)];
			int specie2_count = state[CUDA_GET_SPI(specie_index2, subvolume_index, subvolumes_count)];
#ifdef DEBUG
			printf("hit bi_diff. Specie_index1 = %d, specie_index2 = %d,\n", specie_index, specie_index2);
			printf("    specie1_count = %d, specie2_count = %d\n", specie1_count, specie2_count);
			printf("----------   end react_rate( ) ----------\n\n");
#endif
			return specie1_count * specie2_count * reaction_rate_constants[reaction_number];
		}
	}

	// uni reaction type
	int specie_count = state[CUDA_GET_SPI(specie_index, subvolume_index, subvolumes_count)];
#ifdef DEBUG
	printf("hit uni. Specie_index = %d, specie_count = %d, ", specie_index, specie_count);
	printf("----------   end react_rate( ) ---------- \n\n");
#endif
	return specie_count * reaction_rate_constants[reaction_number];
}

__device__ float * react_rates(int * reactants, int reactions_count, int * state, int subvolumes_count,
		int species_count, int subvolume_index, float * reaction_rate_constants)
{
	__shared__ extern float react_rates_array[];    // we need extern because the size is not a compile-time constant
													// we'll need to allocate during the kernel invocation
													// TODO: rethink about this

	for (int i = 0; i < reactions_count; i++) {
		react_rates_array[i] = react_rate(reactants, reactions_count, state, species_count, subvolumes_count,
				subvolume_index, reaction_rate_constants, i);
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


