#include "../include/cuda/init.cuh"

__global__ void init_prng(curandStateMRG32k3a * prngstate)
{
	INDCHECK()

#ifndef TEST
	curand_init(clock64() * sbi, 0, 0, &prngstate[sbi]);
#else // fixed seeds for a deterministic simulation
	curand_init(sbi, 0, 0, &prngstate[sbi]);
#endif

}

__global__ void init_ncount(neigh neigh)
{
	INDCHECK()

	int nc = 0;
	for (int i = 0; i < 6; i++)
		nc += (neigh.index[sbi * 6 + i] != sbi);
	
	neigh.count[sbi] = nc;

}

__device__ int HOR(reactions reactions, int spi)
{
	int max_hor = 0;
	bool is_bi_reaction = false;

	for (int ri = 0; ri < RC; ri++) {
		// if spi is not a reactant of the current reaction, continue
		// with the next one.
		if (reactions.r[GET_COEFF(spi, ri)] == 0)
			continue;

		// sum all the coeff. of the current reaction to compute its
		// order.
		int hor = 0;
		for (int j = 0; j < SPC; j++) {
			int c = reactions.r[GET_COEFF(j, ri)];
			hor += c;
			// check if ri requires 2 molecules of spi
			if (j == spi && c == 2) {
				is_bi_reaction = true;    // TODO: replace with branchless code
			}
		}

		max_hor = max(hor, max_hor);
	}

	if (is_bi_reaction)
		max_hor = 3;

	return max_hor;
}


__global__ void init_hors(int * hors, reactions reactions, int spc)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi != 0)
		return;

	for (int spi = 0; spi < spc; spi++)
		hors[spi] = HOR(reactions, spi);
}


