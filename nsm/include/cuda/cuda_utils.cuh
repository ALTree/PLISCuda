#ifndef CUDA_UTILS_CUH_
#define CUDA_UTILS_CUH_

// get specie_count index from state using specie_index, subvolume_index and subvolume_count
// use like
//     state[CUDA_GET_SPI(sp_index, sb_index, sbc)];
#define GET_SPI(sp_index, sb_index, sbc) ((sp_index) * (sbc) + (sb_index))

// TODO: move in cuda_utils.cuh
template<typename T>
__device__ T sum_fp_array(T * array, int len)
{
	T sum = 0.0;
	for (int i = 0; i < len; i++)
		sum += array[i];

	return sum;
}

#endif /* CUDA_UTILS_CUH_ */
