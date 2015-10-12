#ifndef NSM_CUH_
#define NSM_CUH_

#include <cuda_runtime.h>
#include <stdio.h>

#include "cuda_utils.cuh"

// TODO: decide a parameters order and stick to it

// computes the react_rate of a single reaction (in one subvolume)
__device__ float react_rate(
		int * state,
		int * reactants,
		int sbc, int spc, int rc,
		int sbi, int ri,
		float * rrc
	);


// returns len(reactions) array with the reactions_rates of all the reactions (in one subvolume)

/*
__device__ float * react_rates(    //
		int * reactants,
		int reactions_count,
		int * state,
		int subvolumes_count,
		int species_count,
		int subvolume_index,    //
		float * reaction_rate_constants    //
		);
*/

__device__ float * react_rates(
		int * state,
		int * reactants,
		int sbc, int spc, int rc,
		int sbi,
		float * rrc
	);

// len(species) array with the diffusion_rates of all the species (in one subvolume)
__device__ float * diff_rates(    //
		int * state,              //
		int subvolumes_count,     //
		int species_count,        //
		int subvolume_index,      //
		float * diffusion_rates_constants    //
		);

__device__ void rate_matrix_row(    //
		int * state,              //
		int * reactants,          //
		int subvolumes_count,     //
		int species_count,        //
		int reactions_count,      //
		float * reaction_rate_constants,     //
		float * diffusion_rate_constants,     //
		float * rate_matrix,       //
		int subvolume_index        //
		);

template<typename T> __device__ T sum_fp_array(T * array, int len);

#endif /* NSM_CUH_ */
