#include "../../include/cuda/leap.cuh"

__device__ bool is_critical(int * state, int * reactants, int * products, int sbi, int ri)
{
	bool crit = false;
	for (int i = 0; i < SPC; i++) {
		crit = crit || ((reactants[GET_COEFF(i, ri)] - products[GET_COEFF(i, ri)]) * NC > state[GET_SPI(i, sbi)]);
	}

	return crit;
}

__device__ float compute_g(int * state, int * reactants, int sbi, int spi)
{

	int hor = HOR(reactants, spi);

#ifdef DEBUG
	if (hor < 1 || hor > 2) {
		printf(">>>>> ARGH! @ compute_g(%d, %d): hor = %d\n", sbi, spi, hor);
	}
#endif

	switch (hor) {
	case 1:
		return 1;
	case 2:
		return 2;
	case 3:
		int x = state[GET_SPI(spi, sbi)];
		if (x == 1) { // TODO: is 1.0 / +Inf == 0? can we use this to avoid the check?
			return 2.0;
		}
		return 2.0 + 1.0 / (state[GET_SPI(spi, sbi)]);
	}
}

__device__ int HOR(int * reactants, int spi)
{
	max_hor = 0;
	bool is_bi_reaction = false;

	for (int i = 0; i < RC; i++) {

		// if spi is not a reactant of the current
		// reaction, continue with the next one.
		if (reactants[GET_COEFF(spi, i)] == 0) {
			continue;
		}

		// sum all the coeff. of the current
		// reaction to compute its order
		hor = 0;
		for (int j = 0; j < SP; j++) {
			int c = reactants[GET_COEFF(j, i)]
			hor += c;
			// check if ri requires 2 molecules of spi
			if (j == spi && c == 2) {
				is_bi_reaction = true;
			}

		}

		max_hor = max(hor, max_hor);
	}

	if (is_bi_reaction) {
		max_hor = 3;
	}

	return max_hor;

}

