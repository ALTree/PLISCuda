#ifndef LEAP_CUH_
#define LEAP_CUH_

#include <cuda_runtime.h>

#include "constants.cuh"
#include "cuda_utils.cuh"

// returns true iff reaction ri is critical in subvolume sbi
__device__ bool is_critical(int * state, int * reactants, int * products, int sbi, int ri);

#endif /* LEAP_CUH_ */
