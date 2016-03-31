#include "../include/cuda/log.cuh"

#include <stdio.h>

__global__ void log_data(int * state, unsigned int * sbv_to_log, int sbv_to_log_len, bool * spc_to_log,
		int spc_to_log_len, int * log_data)
{
	unsigned int ti = blockIdx.x * blockDim.x + threadIdx.x;
	if (ti >= sbv_to_log_len)
		return;

	int spc_index = 0;
	for (int spi = 0; spi < SPC; spi++) {
		if (!spc_to_log[spi])
			continue;
		log_data[spc_index * spc_to_log_len + ti] = state[GET_SPI(spi, sbv_to_log[ti])];
		spc_index++;
	}

}
