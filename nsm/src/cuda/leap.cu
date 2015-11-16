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

	int x = 0;
	switch (hor) {
	case 1:
		return 1;
	case 2:
		return 2;
	case 3:
		x = state[GET_SPI(spi, sbi)];
		if (x == 1) { // TODO: is 1.0 / +Inf == 0? can we use this to avoid the check?
			return 2.0;
		}
		return 2.0 + 1.0 / (x - 1);
	default:
		return 0; // nope
	}
}

__device__ int HOR(int * reactants, int spi)
{
	int max_hor = 0;
	bool is_bi_reaction = false;

	for (int i = 0; i < RC; i++) {

		// if spi is not a reactant of the current
		// reaction, continue with the next one.
		if (reactants[GET_COEFF(spi, i)] == 0) {
			continue;
		}

		// sum all the coeff. of the current
		// reaction to compute its order.
		int hor = 0;
		for (int j = 0; j < SPC; j++) {
			int c = reactants[GET_COEFF(j, i)];
			hor += c;
			// check if ri requires 2 molecules of spi
			if (j == spi && c == 2) {
				is_bi_reaction = true; // TODO: replace with branchless code
			}
		}

		max_hor = max(hor, max_hor);
	}

	if (is_bi_reaction) {
		max_hor = 3;
	}

	return max_hor;

}

__device__ float compute_mu(int * state, int * reactants, int * products, int sbi, int spi, float * react_rates_array)
{
	float mu = 0.0;

	for (int i = 0; i < RC; i++) {

		// when computing mu we only sum over non-critical reactions
		if(is_critical(state, reactants, products, sbi, i)) {
			continue;
		}

		// mu is the sum of (change_vector) * (reaction_rate) over
		// non-critical reactions.
		mu += (products[GET_SPI(spi, sbi)] - reactants[GET_SPI(spi, sbi)]) * react_rates_array[GET_RR(i, sbi)];
	}

	return mu;
}

__device__ float compute_sigma2(int * state, int * reactants, int * products, int sbi, int spi, float * react_rates_array)
{
	float sigma2 = 0.0;

	for (int i = 0; i < RC; i++) {

		// when computing sigma2 we only sum over non-critical reactions
		if(is_critical(state, reactants, products, sbi, i)) {
			continue;
		}

		// sigma2 is the sum of (change_vector)^2 * (reaction_rate) over
		// non-critical reactions.
		int v = products[GET_SPI(spi, sbi)] - reactants[GET_SPI(spi, sbi)];
		sigma2 += (v*v) * react_rates_array[GET_RR(i, sbi)];
	}

	return sigma2;
}


