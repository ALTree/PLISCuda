#ifndef TAU_CUH_
#define TAU_CUH_

#include <cuda_runtime.h>
#include <curand_kernel.h>

#include <math.h>
#include <stdio.h>

#include "constants.cuh"

// Returns true iff reaction ri is critical in subvolume sbi.
__device__ bool is_critical_reaction(int * state, reactions reactions, int sbi, int ri);

// Returns true iff specie spi is critical in subvolume sbi.
__device__ bool is_critical_diffusion(int * state, int sbi, int spi);

// Returns g (as defined in Cao, Gillespie, Petzold - Efficient step
// size selection for the tau-leaping simulation method, J chem Phys
// 124, 044109, page 6) for specie spi in reaction ri inside subvolume
// sbi.
__device__ float compute_g(int * state, int * hors, int sbi, int spi);

// Compute the tau time (as defined in Cao, Gillespie, Petzold -
// Efficient step size selection for the tau-leaping simulation
// method, J chem Phys 124, 044109, page 7, formula 33) for a single
// specie, in subvolume sbi.
// 
// Diffusion events are taken into account (so it's not really [Tao06],
// it's the modified version found in [Harris10]
__device__ float compute_tau_sp(int * state, reactions reactions, int * hors, 
								bool crit_r[MAXREACTIONS], unsigned int * topology,
								int sbi, int spi, rates rates);
	
// Compute the subvolume tau time (i.e. the min of the tau_sp over all
// the species), in subvolume sbi.  The min has to be taken over the
// reactant species NOT involved in critical reactions.  If every
// reaction is critical, returns +Inf.
__device__ float compute_tau_ncr(int * state, reactions reactions, 
								 int * hors, bool crit_r[MAXREACTIONS], unsigned int * topology, 
								 int sbi, rates rates);

// Compute the subvolume tau_cr time (i.e. the tau for the critical
// reactions).  Returns +Inf if every reaction is non-critical.
__device__ float compute_tau_cr(int * state, bool crit_r[MAXREACTIONS],
								int sbi, rates rates, curandStateMRG32k3a * s);

// Fill the tau array with taus computed as [Cao06]
__global__ void compute_taus(int * state, reactions reactions, int * hors, 
							 unsigned int * topology, rates rates, float * tau, 
							 float min_tau, char * leap, curandStateMRG32k3a * s);


#endif /* TAU_CUH_ */
