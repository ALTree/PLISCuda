#ifndef LEAP_CUH_
#define LEAP_CUH_

#include <cuda_runtime.h>
#include <curand_kernel.h>

#include <math.h>
#include <stdio.h>

#include "constants.cuh"
#include "cuda_utils.cuh"


// returns true iff reaction ri is critical in subvolume sbi
__device__ bool is_critical_reaction(int * state, int * reactants, int * products, int sbi, int ri);

// returns true iff specie spi is critical in subvolume sbi
__device__ bool is_critical_diffusion(int * state, int sbi, int spi);

// returns g (as defined in Cao, Gillespie, Petzold - Efficient step size selection
// for the tau-leaping simulation method, J chem Phys 124, 044109, page 6) for specie
// spi in reaction ri inside subvolume sbi
__device__ float compute_g(int * state, int * reactants, int * hors, int sbi, int spi);

// returns HOR(spi). Well, actually it returns
// 1 for HOR(spi) = 1
// 2 for HOR(spi) = 2 and is a "1 1" reaction
// 3 for HOR(spi) = 3 and is a "2" reaction
__device__ int HOR(int * reactants, int spi);

// compute the tau time (as defined in Cao, Gillespie, Petzold - Efficient step size selection
// for the tau-leaping simulation method, J chem Phys 124, 044109, page 7, formula 33) for
// a single specie, in subvolume sbi.
// Diffusion events are taken into account(so it's not really [Tao06], it's che modified version
// found in [Harris10]
__device__ float compute_tau_sp(int * state, int * reactants, int * products, int * hors, 
								bool crit_r[MAXREACTIONS], unsigned int * topology,
								int sbi, int spi, rates rates);
	
// compute the subvolume tau time (i.e. the min of the tau_sp over all the species), in subvolume sbi.
// The min has to be taken over the reactant species NOT involved in critical reactions.
// If every reaction is critical, returns +Inf.
__device__ float compute_tau_ncr(int * state, int * reactants, int * products, 
								 int * hors, bool crit_r[MAXREACTIONS], unsigned int * topology, 
								 int sbi, rates rates);

// compute the subvolume tau_cr time (i.e. the tau for the critical reactions).
// Returns +Inf if every reaction is non-critical.
__device__ float compute_tau_cr(int * state, int * reactants, int * products, bool crit_r[MAXREACTIONS],
								int sbi, rates rates, curandStateMRG32k3a * s);

// Fill the tau array with taus computed as [Cao06]
__global__ void fill_tau_array_leap(int * state, int * reactants, int * products, int * hors, 
									unsigned int * topology, rates rates, float * tau, 
									float min_tau, char * leap, curandStateMRG32k3a * s);

// Performs a single leap step.
// Returns immediately if:
//     - the current thread is associated with a sbv not marked as leap
//
// The kernel fire all non-critical reactions and diffuse all the
// species, and then it fires ONE critical reaction (if it is enabled
// in the current state, after the leaping. If not, it does not fire
// any reaction).  It does NOT update tau or the current time.
__global__ void leap_step(int * state, int * reactants, int * products, unsigned int * topology,
						  rates rates, float min_tau, float * current_time, 
						  char * leap, curandStateMRG32k3a * prngstate);

__global__ void check_state(int * state, bool * revert);

__global__ void initialize_hors_array(int * hors, int * reactants, int spc);

#endif /* LEAP_CUH_ */
