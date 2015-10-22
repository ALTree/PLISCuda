#ifndef NSM_CUH_
#define NSM_CUH_

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>

// fill the whole tau array with random times
__global__ void fill_tau_array(float * tau, int sbc);

// returns the index of a random reaction to fire in the associated subvolume.
// The chance that we choose reaction R is given by the react rate of R over
// the sum of the reaction rates of all the reactions.
__device__ int choose_rand_reaction(int sbc, int rc, float * rate_matrix, float * react_rates_array, float rand);

// returns the index of a random specie to diffuse in the associated subvolume.
// The chance that we choose specie S is given by the diffusion rate of S over
// the sum of the diffusion rates of all the species.
__device__ int choose_rand_specie(int sbc, int spc, float * rate_matrix, float * diff_rates_array, float rand);

void h_fill_tau_array(thrust::device_vector<float> &tau);
int h_get_min_tau(thrust::device_vector<float> &tau);


#endif /* NSM_CUH_ */
