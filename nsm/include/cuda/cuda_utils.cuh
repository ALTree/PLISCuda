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
#define GET_COEFF(spi, ri) ((spi)* (RC) + (ri))

// returns the sum of a floating point array of length len
template<typename T>
__device__ T sum_fp_array(T * array, int len)
{
	T sum = 0.0;
	for (int i = 0; i < len; i++)
		sum += array[i];

	return sum;
}


#endif /* CUDA_UTILS_CUH_ */
