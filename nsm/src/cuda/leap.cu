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
		if (x == 1) {    // TODO: is 1.0 / +Inf == 0? can we use this to avoid the check?
			return 2.0;
		}
		return 2.0 + 1.0 / (x - 1);
	default:
		// HOR(spi) == 0 if spi does not appear as reactant in any reaction.
		// Return +Inf so that when we divide by g we get 0 and the procedure
		// takes max(0, 1) = 1 as g.
		return INFINITY;
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
				is_bi_reaction = true;    // TODO: replace with branchless code
			}
		}

		max_hor = max(hor, max_hor);
	}

	if (is_bi_reaction) {
		max_hor = 3;
	}

	return max_hor;

}

__device__ float compute_mu(int * state, int * reactants, int * products, unsigned int * topology, int sbi, int spi,
		float * react_rates_array, float * diff_rates_array)
{
	float mu = 0.0;

	// sum propensities for the reactions
	for (int i = 0; i < RC; i++) {

		// when computing mu we only sum over non-critical reactions
		if (is_critical(state, reactants, products, sbi, i)) {
			continue;
		}

		// mu is the sum of (change_vector) * (reaction_rate) over
		// non-critical reactions.
		int v = products[GET_COEFF(spi, i)] - reactants[GET_COEFF(spi, i)];
		mu += v * react_rates_array[GET_RR(i, sbi)];
	}

	// add propensities of outgoing diffusions for specie spi.
	// We should sum the diffusion propensities over all the
	// neighbours, but diff_rates_array already has the
	// overall diffusion propensity.
	mu += diff_rates_array[GET_DR(spi, sbi)];

	// add propensities of incoming diffusions for specie spi.
	for (int i = 0; i < 6; i++) {    // loop over the neighbours
		unsigned int ni = topology[sbi * 6 + i];    // neighbour index

		// first we need to compute how many neighbours ni has
		int nni = 0;
		for (int j = 0; j < 6; j++) {
			if (topology[ni * 6 + j] != ni) {
				nni++;
			}
		}

		// now we subtract from mu the propensity of specie spi in
		// subvolume ni divided by nni (i.e. we sum a negative value).
		mu -= (diff_rates_array[GET_DR(spi, ni)]) / nni;

	}

	return mu;
}

__device__ float compute_sigma2(int * state, int * reactants, int * products, unsigned int * topology, int sbi, int spi,
		float * react_rates_array, float * diff_rates_array)
{
	float sigma2 = 0.0;

	for (int i = 0; i < RC; i++) {

		// when computing sigma2 we only sum over non-critical reactions
		if (is_critical(state, reactants, products, sbi, i)) {
			continue;
		}

		// sigma2 is the sum of (change_vector)^2 * (reaction_rate) over
		// non-critical reactions.
		int v = products[GET_COEFF(spi, i)] - reactants[GET_COEFF(spi, i)];
		sigma2 += (v * v) * react_rates_array[GET_RR(i, sbi)];
	}

	// add propensities of outgoing diffusions for specie spi.
	// We should sum the diffusion propensities over all the
	// neighbours, but diff_rates_array already has the
	// overall diffusion propensity.
	sigma2 += diff_rates_array[GET_DR(spi, sbi)];

	// add propensities of incoming diffusions for specie spi.
	for (int i = 0; i < 6; i++) {    // loop over the neighbours
		unsigned int ni = topology[sbi * 6 + i];    // neighbour index

		// first we need to compute how many neighbours ni has
		int nni = 0;
		for (int j = 0; j < 6; j++) {
			if (topology[ni * 6 + j] != ni) {
				nni++;
			}
		}

		// now we add to mu the propensity of specie spi in
		// subvolume ni divided by nni. No need to square since
		// the coeff. is always -1, just sum 1.
		sigma2 += (diff_rates_array[GET_DR(spi, ni)]) / max(nni, 1);    // TODO: fix?

	}

	return sigma2;
}

__device__ float compute_tau_sp(int * state, int * reactants, int * products, unsigned int * topology, int sbi, int spi,
		float * react_rates_array, float * diff_rates_array)
{
	float g = compute_g(state, reactants, sbi, spi);
	int x = state[GET_SPI(spi, sbi)];

	float mu = compute_mu(state, reactants, products, topology, sbi, spi, react_rates_array, diff_rates_array);
	float sigma2 = compute_sigma2(state, reactants, products, topology, sbi, spi, react_rates_array, diff_rates_array);

	float m = max(EPSILON * x / g, 1.0f);
	float t1 = m / abs(mu);
	float t2 = (m * m) / (sigma2);

	return min(t1, t2);
}

__device__ float compute_tau(int * state, int * reactants, int * products, unsigned int * topology, int sbi,
		float * react_rates_array, float * diff_rates_array)
{
	float min_tau = INFINITY;

	for (int spi = 0; spi < SPC; spi++) {

		// First of all we need to check if the specie spi is involved
		// as reactant in a critical reaction. If it is, skip it.
		bool skip = false;
		for (int ri = 0; ri < RC; ri++) {    // iterate over reactions
			if (is_critical(state, reactants, products, sbi, ri)) {    // if it's critical
				// skip if the specie spi is involved in the reaction
				skip = skip || (reactants[GET_COEFF(spi, ri)] > 0);
			}
		}

		if (skip) {
			continue;
		}
		// spi is not involved in critical reactions.

		float tau = compute_tau_sp(state, reactants, products, topology, sbi, spi, react_rates_array, diff_rates_array);

		min_tau = min(min_tau, tau);
	}

	return min_tau;
}

__global__ void fill_tau_array_leap(int * state, int * reactants, int * products, unsigned int * topology,
		float * rate_matrix, float * react_rates_array, float * diff_rates_array, float * tau, bool * leap)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	float tauv = compute_tau(state, reactants, products, topology, sbi, react_rates_array, diff_rates_array);
	tau[sbi] = tauv;
	leap[sbi] = tauv >= 1.0 / rate_matrix[GET_RATE(2, sbi)];
}

__global__ void leap_step(int * state, int * reactants, int * products, unsigned int * topology,
		float * react_rates_array, float * diff_rates_array, float * tau, bool * leap, curandStateMRG32k3a * prngstate)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC || !leap[sbi])
		return;

	// fire all the reaction events
	for (int ri = 0; ri < RC; ri++) {

		unsigned int k;    // how many times it fires
		k = curand_poisson(&prngstate[sbi], tau[sbi] * react_rates_array[GET_RR(ri, sbi)]);

		// update state
		for (int spi = 0; spi < SPC; spi++) {
			state[GET_SPI(spi, sbi)] += k * (products[GET_COEFF(spi, ri)] - reactants[GET_COEFF(spi, ri)]);
		}

	}

	__syncthreads();

}

