#ifndef NSM_CUH_
#define NSM_CUH_

#include <cuda_runtime.h>
#include <stdio.h>

#include "cuda_utils.cuh"

// computes the react_rate of a single reaction
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

__device__ float * react_rates(    //
		int * reactants,    //
		int reactions_count,    //
		int * state,    //
		int subvolumes_count,    //
		int species_count,    //
		int subvolume_index,    //
		float * reaction_rate_constants    //
		);

__device__ float sum_react_rates(float * react_rates, int reactions_count);

__device__ float * diff_rates(    //
		int * state,              //
		int subvolumes_count,     //
		int species_count,        //
		int subvolume_index,      //
		float * diffusion_rates_constants  //
);

__device__ float sum_diff_rates(float * diff_rates, int species_count);

#endif /* NSM_CUH_ */
