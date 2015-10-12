#ifndef CUDA_UTILS_CUH_
#define CUDA_UTILS_CUH_

// get specie_count index from specie_index, subvolume_index and subvolume_count
// use like
//     state[CUDA_GET_SPI(sp_index, sb_index, sbc)];
#define CUDA_GET_SPI(sp_index, sb_index, sbc) ((sp_index) * (sbc) + (sb_index))

#endif /* CUDA_UTILS_CUH_ */
