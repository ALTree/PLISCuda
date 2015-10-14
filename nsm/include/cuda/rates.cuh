#ifndef NSM_CUH_
#define NSM_CUH_

#include <cuda_runtime.h>
#include <stdio.h>

#include "cuda_utils.cuh"

// computes the react_rate of a single reaction (in one subvolume)
__device__ float react_rate(int * state, int * reactants, int sbc, int spc, int rc, int sbi, int ri, float * rrc);

// computes the reactions rates of all the reactions (in all the subvolumes)
// writes the result (rc * sbc float values) at the address in global memory pointed by react_rates_array
__device__ void react_rates(int * state, int * reactants, int sbc, int spc, int rc, float * rrc,
		float * react_rates_array);

// computes the diffusion rates of all the reactions (in all the subvolumes)
// writes the result (rc * spc float values) at the address in global memory pointed by diff_rates_array
__device__ void diff_rates(int * state, int sbc, int spc, float * drc, float * diff_rates_array);

// updated a single row of the rate matrix
__device__ void rate_matrix_row(int * state, int * reactants, int sbc, int spc, int rc, int sbi, float * rate_matrix,
		float * rrc, float * drc);

// returns the sum of a floating point array of length len
template<typename T> __device__ T sum_fp_array(T * array, int len);

#endif /* NSM_CUH_ */
