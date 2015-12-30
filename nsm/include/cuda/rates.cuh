#ifndef RATES_CUH_
#define RATES_CUH_

#include <cuda_runtime.h>
#include <stdio.h>

#include "cuda_utils.cuh"
#include "constants.cuh"

// computes the react_rate of a single reaction (in one subvolume)
__device__ float react_rate(int * state, int * reactants, int sbi, int ri, float * rrc);

// computes the reactions rates of all the reactions (in all the subvolumes)
// writes the result (RC * SBC float values) at the address in global memory pointed by react_rates_array
__device__ void react_rates(int * state, int * reactants, float * rrc, float * react_rates_array);

// computes the diffusion rates of all the reactions (in all the subvolumes)
// writes the result (RC * SPC float values) at the address in global memory pointed by diff_rates_array
__device__ void diff_rates(int * state, float * drc, float * diff_rates_array);

// updates the whole rate matrix
__device__ void update_rate_matrix(unsigned int * topology, float * rate_matrix, float * react_rates_array,
		float * diff_rates_array);

// call react_rates, diff_rates and update_rate_matrix in the correct way
__global__ void compute_rates(int * state, int * reactants, unsigned int * topology, float * rate_matrix, float * rrc,
		float * drc, int * d_subv_consts, float * react_rates_array, float * diff_rates_array);

#endif /* RATES_CUH_ */
