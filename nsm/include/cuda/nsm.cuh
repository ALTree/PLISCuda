#ifndef NSM_CUH_
#define NSM_CUH_

#include <cuda_runtime.h>
#include <stdio.h>

#include "cuda_utils.cuh"

// TODO: decide a parameters order and stick to it

// computes the react_rate of a single reaction (in one subvolume)
__device__ float react_rate(    //
		int * reactants,          // reactants array
		int reactions_count,       //
		int * state,              // state array
		int species_count,         // the number of species
		int subvolumes_count,       // the number of subvolumes
		int subvolume_index,      // the subvolume we are operating in
		float * reaction_rate_constants,     // reaction rate constants array
		int reaction_number       // the reaction we are processing
		);

// len(reactions) array with the reactions_rates of all the reactions (in one subvolume)
__device__ float * react_rates(    //
		int * reactants,    //
		int reactions_count,    //
		int * state,    //
		int subvolumes_count,    //
		int species_count,    //
		int subvolume_index,    //
		float * reaction_rate_constants    //
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
