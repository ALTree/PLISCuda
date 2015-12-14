#include "../../include/cuda/leap.cuh"

__device__ bool is_critical_reaction(int * state, int * reactants, int * products, int sbi, int ri)
{
	bool crit = false;
	for (int spi = 0; spi < SPC; spi++) {
		crit = crit || ((reactants[GET_COEFF(spi, ri)] - products[GET_COEFF(spi, ri)]) * NC > state[GET_SPI(spi, sbi)]);
	}

	return crit;
}

__device__ bool is_critical_diffusion(int * state, int sbi, int spi)
{
	return state[GET_SPI(spi, sbi)] < NC;
}

__device__ float compute_g(int * state, int * reactants, int sbi, int spi)
{

	int hor = HOR(reactants, spi);

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
		if (is_critical_reaction(state, reactants, products, sbi, i)) {
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
		if (is_critical_reaction(state, reactants, products, sbi, i)) {
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

__device__ float compute_tau_ncr(int * state, int * reactants, int * products, unsigned int * topology, int sbi,
		float * react_rates_array, float * diff_rates_array)
{
	float min_tau = INFINITY;

	for (int spi = 0; spi < SPC; spi++) {

		// First of all we need to check if the specie spi is involved
		// in a critical event. If it is, skip it.
		bool skip = false;

		// check for critical reaction events
		for (int ri = 0; ri < RC; ri++) {    // iterate over reactions
			if (is_critical_reaction(state, reactants, products, sbi, ri)) {    // if it's critical
				// skip if the specie spi is involved in the reaction
				skip = skip || (reactants[GET_COEFF(spi, ri)] > 0);
			}
		}

		// check for critical diffusion events
		skip = skip || is_critical_diffusion(state, sbi, spi);

		if (skip) {
			continue;
		}
		// spi is not involved in any critical event.

		float tau = compute_tau_sp(state, reactants, products, topology, sbi, spi, react_rates_array, diff_rates_array);

		min_tau = min(min_tau, tau);
	}

	return min_tau;
}

__device__ float compute_tau_cr(int * state, int * reactants, int * products, int sbi, float * react_rates_array,
		float * diff_rates_array, curandStateMRG32k3a * s)
{
	float react_rates_sum_cr = 0.0;    // sum of the react rates of critical reactions
	float diff_rates_sum_cr = 0.0;     // sum of diffusion rates of critical diffusion events

	for (int ri = 0; ri < RC; ri++) {
		react_rates_sum_cr += (react_rates_array[GET_RR(ri, sbi)]
				* is_critical_reaction(state, reactants, products, sbi, ri));
	}

	for (int spi = 0; spi < SPC; spi++) {
		diff_rates_sum_cr += (diff_rates_array[GET_DR(spi, sbi)] * is_critical_diffusion(state, sbi, spi));
	}

	if (react_rates_sum_cr == 0.0 && diff_rates_sum_cr == 0.0)
		return INFINITY;

	float rand = curand_uniform(&s[sbi]);
	return -logf(rand) / (react_rates_sum_cr + diff_rates_sum_cr);
}

__global__ void fill_tau_array_leap(int * state, int * reactants, int * products, unsigned int * topology,
		float * rate_matrix, float * react_rates_array, float * diff_rates_array, float * tau, char * leap,
		curandStateMRG32k3a * s)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	float tau_ncr = compute_tau_ncr(state, reactants, products, topology, sbi, react_rates_array, diff_rates_array);
	float tau_cr = compute_tau_cr(state, reactants, products, sbi, react_rates_array, diff_rates_array, s);

	// if tau_ncr is too small, we can't leap in this subvolume.
	// Also Prevent leap = true if tau_ncr is +Inf
	bool leap_here = true;
	if (isinf(tau_ncr) || (tau_ncr < 2.0 / rate_matrix[GET_RATE(2, sbi)])) {
		leap[sbi] = SSA;
		leap_here = false;
	}

	if (tau_ncr < tau_cr) {
		// no critical event will happen, we'll leap with
		// all the non-critical events
		tau[sbi] = tau_ncr;
		if (leap_here)
			leap[sbi] = LEAP_NOCR;
	} else {
		// a single critical event will happen, all the
		// non-critical events will leap with tau = tau_cr
		tau[sbi] = tau_cr;
		if (leap_here)
			leap[sbi] = LEAP_CR;
	}

}

__global__ void leap_step(int * state, int * reactants, int * products, float * rate_matrix, unsigned int * topology,
		float * react_rates_array, float * diff_rates_array, float * rrc, float * drc, float min_tau,
		float * current_time, char * leap, curandStateMRG32k3a * prngstate)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC || leap[sbi] == SSA)
		return;

	// count neighbours of the current subvolume. We'll need the value later.
	// TODO: remove when issue #26 is fixed
	int neigh_count = 0;
	for (int i = 0; i < 6; i++)
		neigh_count += (topology[sbi * 6 + i] != sbi);

	// fire all the non-critical reaction events
	for (int ri = 0; ri < RC; ri++) {

		if (is_critical_reaction(state, reactants, products, sbi, ri)) {
			continue;
		}

		unsigned int k;    // how many times it fires
		k = curand_poisson(&prngstate[sbi], min_tau * react_rates_array[GET_RR(ri, sbi)]);

		if (k > 0)
			printf("(%f) [subv %d] fire reaction %d for %d times\n", *current_time, sbi, ri, k);

		// update state
		// TODO: needs to be atomic? I suspect so..
		// Maybe not if we use __synchthreads( ) before the diffusion events.
		for (int spi = 0; spi < SPC; spi++) {
			int new_state = k * (products[GET_COEFF(spi, ri)] - reactants[GET_COEFF(spi, ri)]);
			atomicAdd(&state[GET_SPI(spi, sbi)], new_state);
		}

	}

	__syncthreads();

	// fire non-critical outgoing diffusion events
	for (int spi = 0; spi < SPC; spi++) {

		if (is_critical_diffusion(state, sbi, spi)) {
			continue;
		}

		// diffuse to each neighbour
		for (unsigned int ngb = 0; ngb < neigh_count; ngb++) {
			unsigned int k = curand_poisson(&prngstate[sbi], min_tau * diff_rates_array[GET_DR(spi, sbi)]);
			atomicAdd(&state[GET_SPI(spi, topology[sbi*6 + ngb])], k);
			atomicSub(&state[GET_SPI(spi, sbi)], k);
			if (k > 0)
				printf("(%f) [subv %d] diffuse %d molecules of specie %d to subv %d \n", *current_time, sbi, k, spi,
						topology[sbi * 6 + ngb]);
		}

	}

	if (leap[sbi] != LEAP_CR)
		return;

	// Problem: in the following we use the old react_rates_array (we'll
	// update it later, in another kernel call), but it can happen that a
	// reaction has rate > 0 even if we can't fire it (for example because
	// we had 2 molecules at the start of this step, but the leap diffusion
	// phase we just executed removed one of them).
	// To avoid negative population, we choose the random reaction as usual,
	// but we only fire it if we have enough molecules in the subvolume.
	// Note that it's possibile that something has entered from a neighbour
	// during the leap phase, so checking the state is a good pragmatic
	// way to ensure that we'll fire everytime is possibile.

	__syncthreads();

	// fire a single critical event
	float rand = curand_uniform(&prngstate[sbi]);

	// first we have to choose if we'll fire a reaction or we'll diffuse a molecule

	// sum the reaction rates of critical reactions
	float rr_sum = 0.0;
	for (int ri = 0; ri < RC; ri++)
		rr_sum += react_rates_array[GET_RR(ri, sbi)] * is_critical_reaction(state, reactants, products, sbi, ri);

	// sum the diffusion rates of critical diffusion events
	float dr_sum = 0.0;
	for (int spi = 0; spi < SPC; spi++)
		dr_sum += diff_rates_array[GET_DR(spi, sbi)] * is_critical_diffusion(state, sbi, spi);

	// if the sum is zero we can't fire or diffuse anything
	if (rr_sum + dr_sum == 0.0) {
		return;
	}

	if (rand < rr_sum / (rr_sum + dr_sum)) {    // reaction

		float scaled_sum = rr_sum * rand;
		float partial_sum = 0;

		int ric = 0;
		while (partial_sum <= scaled_sum) {
			partial_sum += react_rates_array[GET_RR(ric, sbi)]
					* is_critical_reaction(state, reactants, products, sbi, ric);
			ric++;
		}
		// We'll fire the ric-nth critical reactions.

		int ri;
		for (ri = 0; ri < RC && ric > 0; ri++) {
			if (is_critical_reaction(state, reactants, products, sbi, ri)) {
				ric--;
			}
		}
		ri = ri - 1;

		// Check if the current state lets us fire reaction ri
		// (see comment above).
		bool fire = true;
		for (int spi = 0; spi < SPC; spi++) {
			fire = fire && (state[GET_SPI(spi, sbi)] >= reactants[GET_COEFF(spi, ri)]);
		}

		if (!fire)
			return;

		for (int spi = 0; spi < SPC; spi++) {    // TODO: atomic add?
			state[GET_SPI(spi, sbi)] += (products[GET_COEFF(spi, ri)] - reactants[GET_COEFF(spi, ri)]);
		}

		printf("(%f) [subv %d] fire reaction %d (critical)\n", *current_time, sbi, ri);

	} else {    // diffusion

		float scaled_sum = dr_sum * rand;    // no need to scale down by neigh_count (it's the raw sum)
		float partial_sum = 0;

		int spic = 0;
		while (partial_sum <= scaled_sum) {
			partial_sum += diff_rates_array[GET_DR(spic, sbi)] * is_critical_diffusion(state, sbi, spic);
			spic++;
		}
		// We'll diffuse the spic-nth critical specie.

		int spi;
		for (spi = 0; spi < SPC && spic > 0; spi++) {
			if (is_critical_diffusion(state, sbi, spi)) {
				spic--;
			}
		}
		spi = spi - 1;

		// Check if the current state lets us diffuse specie spi
		// (see comment above).
		bool fire = state[GET_SPI(spi, sbi)] > 0;

		if (!fire)
			return;

		// choose a random destination
		// TODO: we need to re-use the rand we already have.
		int rdi;
		do {
			rdi = (int) (curand_uniform(&prngstate[sbi]) * 6);
		} while (rdi > 5);

		// get index of neighbour #rdi (overwrite rdi, whatever)
		rdi = topology[sbi * 6 + rdi];

		// If rdi == sbi (i.e. diffuse to myself) don't do anything
		// TODO: atomic?
		if (rdi != sbi) {
			atomicAdd(&state[GET_SPI(spi, sbi)], -1);
			atomicAdd(&state[GET_SPI(spi, rdi)], 1);
		}

		printf("(%f) [subv %d] diffuse specie %d to %d (critical)\n", *current_time, sbi, spi, rdi);
	}

}

__global__ void check_state(int * state, bool * revert)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	bool _revert = false;
	for (int spi = 0; spi < SPC; spi++)
		_revert = _revert || (state[GET_SPI(spi, sbi)] < 0);

	revert[sbi] = _revert;
}

