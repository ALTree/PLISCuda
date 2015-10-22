#ifndef NSM_CUH_
#define NSM_CUH_

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>

#include "cuda_utils.cuh"
#include "rates.cuh"

// fill the whole tau array with random times
__global__ void fill_tau_array(float * tau, float * rate_matrix, int sbc);

// returns the index of a random reaction to fire in the associated subvolume.
// The chance that we choose reaction R is given by the react rate of R over
// the sum of the reaction rates of all the reactions.
__device__ int choose_rand_reaction(int sbc, int rc, float * rate_matrix, float * react_rates_array, float rand);

// returns the index of a random specie to diffuse in the associated subvolume.
// The chance that we choose specie S is given by the diffusion rate of S over
// the sum of the diffusion rates of all the species.
__device__ int choose_rand_specie(int * topology, int sbc, int spc, float * rate_matrix, float * diff_rates_array, float rand);

int h_get_min_tau(thrust::device_vector<float> &tau);

__global__ void nsm_step(int * state, int * reactants, int * products, int * topology, int sbc, int spc, int rc,
		float * rate_matrix, float * rrc, float * drc, float * react_rates_array, float * diff_rates_array, float * tau, int min_sbi);

#endif /* NSM_CUH_ */
