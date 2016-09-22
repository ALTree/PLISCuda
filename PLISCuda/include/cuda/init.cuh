#ifndef INIT_CUH_
#define INIT_CUH_

#include <curand_kernel.h>

#include "constants.cuh"

// Initialize prngstate array (we need a different prng for each thread).
__global__ void init_prng(curandStateMRG32k3a * prngstate);

// Initialize the neigh.count array. neigh.count[sbi] containts the
// number of neighbours the subvolume sbi has.
__global__ void init_ncount(neigh neigh);

// returns HOR(spi). Well, actually it returns
// 1 for HOR(spi) = 1
// 2 for HOR(spi) = 2 and is a "1 1" reaction
// 3 for HOR(spi) = 3 and is a "2" reaction
__device__ int HOR(reactions reactions, int spi);

// Initialize the hors array. Since the HORs values never change, it
// makes sense to cache them.
__global__ void init_hors(int * hors, reactions reactions, int spc);

#endif /* INIT_CUH_ */
