#ifndef LEAP_CUH_
#define LEAP_CUH_

#include <cuda_runtime.h>
#include <curand_kernel.h>

#include <math.h>
#include <stdio.h>

#include "constants.cuh"
#include "tau.cuh"

// Performs a single leap step.
// Returns immediately if:
//     - the current thread is associated with a sbv not marked as leap
//
// The kernel fire all non-critical reactions and diffuse all the
// species, and then it fires ONE critical reaction (if it is enabled
// in the current state, after the leaping. If not, it does not fire
// any reaction).  It does NOT update tau or the current time.
__global__ void leap_step(state state, reactions reactions, unsigned int * topology,
						  rates rates, float min_tau, float * current_time, 
						  char * leap, curandStateMRG32k3a * prngstate);

__global__ void check_state(state state, int * revert);


// returns HOR(spi). Well, actually it returns
// 1 for HOR(spi) = 1
// 2 for HOR(spi) = 2 and is a "1 1" reaction
// 3 for HOR(spi) = 3 and is a "2" reaction
__device__ int HOR(reactions reactions, int spi);

__global__ void initialize_hors_array(int * hors, reactions reactions, int spc);


#endif /* LEAP_CUH_ */
