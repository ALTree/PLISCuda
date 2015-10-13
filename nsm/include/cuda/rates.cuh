#ifndef NSM_CUH_
#define NSM_CUH_

#include <cuda_runtime.h>
#include <stdio.h>

#include "cuda_utils.cuh"

// TODO: decide a parameters order and stick to it

// computes the react_rate of a single reaction (in one subvolume)
__device__ float react_rate(int * state, int * reactants, int sbc, int spc, int rc, int sbi, int ri, float * rrc);

// returns len(reactions) array with the reactions_rates of all the reactions (in one subvolume)
__device__ float * react_rates(int * state, int * reactants, int sbc, int spc, int rc, int sbi, float * rrc);

// len(species) array with the diffusion_rates of all the species (in one subvolume)
__device__ float * diff_rates(int * state, int sbc, int spc, int sbi, float * drc);

// updated a single row of the rate matrix
__device__ void rate_matrix_row(int * state, int * reactants, int sbc, int spc, int rc, int sbi, float * rate_matrix,
		float * rrc, float * drc);

// returns the sum of a floating point array of length len
template<typename T> __device__ T sum_fp_array(T * array, int len);

#endif /* NSM_CUH_ */
