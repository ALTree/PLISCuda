#ifndef RATES_CUH_
#define RATES_CUH_

#include <cuda_runtime.h>
#include <stdio.h>

#include "cuda_utils.cuh"
#include "constants.cuh"

// computes the react_rate of a single reaction (in one subvolume)
__device__ float react_rate(int * state, int * reactants, int sbi, int ri, rates rates);

// computes the reactions rates of all the reactions (in all the subvolumes)
// writes the result (RC * SBC float values) at the address in global memory pointed by react_rates_array
__device__ void react_rates(int * state, int * reactants, rates rates);

// computes the diffusion rates of all the reactions (in all the subvolumes)
// writes the result (RC * SPC float values) at the address in global memory pointed by diff_rates_array
__device__ void diff_rates(int * state, rates rates);

// updates the whole rate matrix
__device__ void update_rate_matrix(unsigned int * topology, rates rates);

// call react_rates, diff_rates and update_rate_matrix in the correct way
__global__ void compute_rates(int * state, int * reactants, unsigned int * topology, 
							  rates rates, int * d_subv_consts);

#endif /* RATES_CUH_ */
