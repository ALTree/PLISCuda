#include "constants.cuh"

#ifndef CUDA_UTILS_CUH_
#define CUDA_UTILS_CUH_

// get specie_count index from state using specie index, subvolume index and subvolume count
// use like
//     state[CUDA_GET_SPI(spi, sbi)];
#define GET_SPI(spi, sbi) ((spi) * (SBC) + (sbi))

// get R or D or R+D index from rate matrix using subvolume index and 0, 1, 2 (resp. R, D, R+D)
// use like
//     rate_matrix[GET_RATE(i, sbi)]
#define GET_RATE(i, sbi) ((i) * (SBC) + (sbi))

// get the stechiometric coefficient for specie spi in reaction ri
// use like
//     reactants[GET_COEFF(spi, ri)]
#define GET_COEFF(spi, ri) ((spi) * (RC) + (ri))

// get the react rate for reaction ri in subvolume sbi
// use like
//     react_rates_array[(GET_RR(ri, sbi)]
#define GET_RR(ri, sbi) ((ri) * (SBC) + (sbi))

// get the diff rate for specie spi in subvolume sbi
// use like
//     diff_rates_array[(GET_DR(spi, sbi)]
#define GET_DR(spi, sbi) ((spi) * (SBC) + (sbi))

#endif /* CUDA_UTILS_CUH_ */
