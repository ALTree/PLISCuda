#ifndef LEAP_CUH_
#define LEAP_CUH_

#include <cuda_runtime.h>

#include <math.h>

#include "constants.cuh"
#include "cuda_utils.cuh"

// returns true iff reaction ri is critical in subvolume sbi
__device__ bool is_critical(int * state, int * reactants, int * products, int sbi, int ri);

// returns g (as defined in Cao, Gillespie, Petzold - Efficient step size selection
// for the tau-leaping simulation method, J chem Phys 124, 044109, page 6) for specie
// spi in reaction ri inside subvolume sbi
__device__ float compute_g(int * state, int * reactants, int sbi, int spi);

// returns HOR(i). Well, actually it returns
// 1 for HOR(i) = 1
// 2 for HOR(i) = 2 and is a "1 1" reaction
// 3 for HOR(i) = 3 and is a "2" reaction
__device__ int HOR(int * reactants, int spi);

// compute mu (as defined in Cao, Gillespie, Petzold - Efficient step size selection
// for the tau-leaping simulation method, J chem Phys 124, 044109, page 7, formula 32a),
// for specie spi in subvolume sbi
__device__ float compute_mu(int * state, int * reactants, int * products, unsigned int * topology, int sbi, int spi,
		float * react_rates_array, float * diff_rates_array);

// compute sigma2 (as defined in Cao, Gillespie, Petzold - Efficient step size selection
// for the tau-leaping simulation method, J chem Phys 124, 044109, page 7, formula 32b),
// for specie spi in subvolume sbi
__device__ float compute_sigma2(int * state, int * reactants, int * products, unsigned int * topology, int sbi, int spi,
		float * react_rates_array, float * diff_rates_array);

// compute the tau time (as defined in Cao, Gillespie, Petzold - Efficient step size selection
// for the tau-leaping simulation method, J chem Phys 124, 044109, page 7, formula 33) for
// a single specie, in subvolume sbi.
// Not to be called if spi is involved in any critical reaction.
__device__ float compute_tau_sp(int * state, int * reactants, int * products, unsigned int * topology, int sbi, int spi,
		float * react_rates_array, float * diff_rates_array);

// compute the subvolume tau time (i.e. the min of the tau_sp over all the species), in subvolume sbi.
// The min has to be taken over the reactant species NOT involved in critical reactions.
// If every reaction is critical, returns +Inf.
__device__ float compute_tau(int * state, int * reactants, int * products, unsigned int * topology, int sbi,
		float * react_rates_array, float * diff_rates_array);

__device__ void fill_tau_array_leap(int * state, int * reactants, int * products, unsigned int * topology,
		float * react_rates_array, float * diff_rates_array, float * tau);

#endif /* LEAP_CUH_ */
