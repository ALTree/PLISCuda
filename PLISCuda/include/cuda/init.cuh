#ifndef INIT_CUH_
#define INIT_CUH_

#include <curand_kernel.h>

#include "constants.cuh"

// initialize prngstate array (one prng for each thread).
__global__ void initialize_prngstate_array(curandStateMRG32k3a * prngstate);

__global__ void init_neigh_count(neigh neigh);

// returns HOR(spi). Well, actually it returns
// 1 for HOR(spi) = 1
// 2 for HOR(spi) = 2 and is a "1 1" reaction
// 3 for HOR(spi) = 3 and is a "2" reaction
__device__ int HOR(reactions reactions, int spi);

__global__ void initialize_hors_array(int * hors, reactions reactions, int spc);

#endif /* INIT_CUH_ */
