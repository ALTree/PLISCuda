#ifndef LEAP_CUH_
#define LEAP_CUH_

#include <cuda_runtime.h>

#include "constants.cuh"
#include "cuda_utils.cuh"

// returns true iff reaction ri is critical in subvolume sbi
__device__ bool is_critical(int * state, int * reactants, int * products, int sbi, int ri);
__device__ float compute_m();

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
__device__ float compute_mu(int * state, int * reactants, int * products, int sbi, int spi, float * react_rates_array);

#endif /* LEAP_CUH_ */
