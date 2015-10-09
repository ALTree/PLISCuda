#ifndef NSM_CUH_
#define NSM_CUH_

#include <cuda_runtime.h>
#include <stdio.h>

// computes the react_rate of a single reaction
__device__ float react_rate(    //
		int * reactants,          // reactants array
		int reactions_count,       //
		int * state,              // state array
		int species_count,         // the number of species
		int subvolumes_count,       // the number of subvolumes
		int subvolume_index,      // the subvolume we are operating in
		double * reaction_rate_constants,     // reaction rate constants array
		int reaction_number       // the reaction we are processing
		);

__device__ float * react_rates(    //
		int * reactants,    //
		int reactions_count,    //
		int * state,    //
		int subvolumes_count,    //
		int species_count,    //
		int subvolume_index,    //
		double * reaction_rate_constants    //
		);



__global__ void test();
void foo();

#endif /* NSM_CUH_ */
