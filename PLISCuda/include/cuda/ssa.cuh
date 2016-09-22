#ifndef NSM_CUH_
#define NSM_CUH_

#include <curand_kernel.h>

#include "rates.cuh"
#include "constants.cuh"

// returns the index of a random reaction to fire in the associated subvolume.
// The chance that we choose reaction R is given by the react rate of R over
// the sum of the reaction rates of all the reactions.
__device__ int choose_rand_reaction(rates rates, float rand);

// returns the index of a random specie to diffuse in the associated subvolume.
// The chance that we choose specie S is given by the diffusion rate of S over
// the sum of the diffusion rates of all the species.
__device__ int choose_rand_specie(neigh neigh, rates rates, float rand);

// Performs a single SSA step.
// Returns immediately if:
//     - the current thread is associated with a sbv marked as leap
//     - the subvolume index is different from min_sbi (the one with minimum tau)
//
// The kernel either fire a reaction or performs a diffusion.
// It updates the state of the neighbours, if necessary.
// It does NOT update tau or the current time.
__global__ void ssa_step(state state, reactions reactions, neigh neigh, rates rates,
						 int min_sbi, float * current_time, char * leap,
						 curandStateMRG32k3a * s);

#endif /* NSM_CUH_ */
