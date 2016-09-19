#ifndef RATES_CUH_
#define RATES_CUH_

#include <cuda_runtime.h>
#include <stdio.h>

#include "constants.cuh"

// computes the react_rate of a single reaction (in one subvolume)
__device__ float react_rate(state state, reactions reactions, int sbi, int ri, rates rates);

// computes the reactions rates of all the reactions (in all the
// subvolumes) writes the result (RC * SBC float values) at the
// address in global memory pointed by react_rates_array
__device__ void react_rates(state state, reactions reactions, rates rates);

// computes the diffusion rates of all the reactions (in all the
// subvolumes) writes the result (RC * SPC float values) at the
// address in global memory pointed by diff_rates_array
__device__ void diff_rates(state state, rates rates);

// updates the whole rate matrix
__device__ void update_rate_matrix(unsigned int * topology, rates rates);

// call react_rates, diff_rates and update_rate_matrix
__global__ void compute_rates(state state, reactions reactions, unsigned int * topology, 
							  rates rates, int * d_subv_consts);

#endif /* RATES_CUH_ */
