#include "../include/cuda/tau.cuh"

__device__ bool is_critical_reaction(state state, reactions reactions, int sbi, int ri)
{
	bool crit = false;
	for (int spi = 0; spi < SPC; spi++) {
		// if state == 0 and the reactions requires specie spi, it's
		// obviously critical
		if(reactions.r[GET_COEFF(spi, ri)] > 0 && state.curr[GET_SPI(spi, sbi)] == 0){
			return true;
		}

		// delta is the net variation in population spi after NC
		// firings of the reaction
		int delta = (reactions.p[GET_COEFF(spi, ri)] - reactions.r[GET_COEFF(spi, ri)]) * NC;
		if(delta >= 0) {
			// the reaction actually *increases* (or leave unchanged)
			// the current specie popolation, so it's obviously not
			// critical
			continue;
		}

		// Now delta < 0 and abs(delta) is the decrease in specie spi
		// population caused by the reaction. 
		// If abs(delta) > population, the reaction is critical.
		crit = crit || (abs(delta) > state.curr[GET_SPI(spi, sbi)]);
	}

	return crit;
}

__device__ bool is_critical_diffusion(state state, int sbi, int spi)
{
	return state.curr[GET_SPI(spi, sbi)] < NC;
}

__device__ float compute_g(state state, reactions reactions, int * hors, int sbi, int spi)
{
	int hor = hors[spi];

	int x = 0;
	switch (hor) {
	case 1:
		return 1;
	case 2:
		return 2;
	case 3:
		x = state.curr[GET_SPI(spi, sbi)];
		if (x == 1) {    // TODO: is 1.0 / +Inf == 0? can we use this to avoid the check?
			return 2.0;
		}
		return 2.0 + 1.0 / (x - 1);
	default:
		// HOR(spi) == 0 if spi does not appear as reactant in any
		// reaction.  Return +Inf so that when we divide by g we get 0
		// and the procedure takes max(0, 1) = 1 as g.
		return INFINITY;
	}
}

__device__ float compute_tau_sp(state state, reactions reactions, int * hors, 
								bool crit_r[MAXREACTIONS], neigh neigh,
								int sbi, int spi, rates rates) 
{
	float g = compute_g(state, reactions, hors, sbi, spi);
	int x = state.curr[GET_SPI(spi, sbi)];

	// compute mu and sigma2

	// compute mu (as defined in Cao, Gillespie, Petzold - Efficient
	// step size selection for the tau-leaping simulation method, J
	// chem Phys 124, 044109, page 7, formula 32a), for specie spi in
	// subvolume sbi
	float mu = 0.0;

	// compute sigma2 (as defined in Cao, Gillespie, Petzold -
	// Efficient step size selection for the tau-leaping simulation
	// method, J chem Phys 124, 044109, page 7, formula 32b), for
	// specie spi in subvolume sbi
	float sigma2 = 0.0;

	// sum propensities for the reactions
	for (int ri = 0; ri < RC; ri++) {

		// when computing mu we only sum over non-critical reactions
		if (crit_r[ri]) {
			continue;
		}

		// mu is the sum of (change_vector) * (reaction_rate) over
		// non-critical reactions.
		// 
		// sigma2 is the sum of (change_vector)² * (reaction_rate)
		// over non-critical reactions.
		int v = reactions.p[GET_COEFF(spi, ri)] - reactions.r[GET_COEFF(spi, ri)];
		mu += v * rates.reaction[GET_RR(ri, sbi)];
		sigma2 += (v * v) * rates.reaction[GET_RR(ri, sbi)];
	}

	if(is_critical_diffusion(state, sbi, spi)) {
		// if spi is critical in this subvolume, don't sum
		// propensities of outgoing diffusions
	} else {
		// Add propensities of outgoing diffusions for specie spi.  We
		// should sum the diffusion propensities over all the
		// neighbours, but diff_rates_array already has the overall
		// diffusion propensity.
		mu += rates.diffusion[GET_DR(spi, sbi)];
		sigma2 += rates.diffusion[GET_DR(spi, sbi)];
	}

	// add propensities of incoming diffusions for specie spi
	for (int i = 0; i < 6; i++) {    // loop over the neighbours
		unsigned int ni = neigh.index[sbi * 6 + i];    // neighbour index
		if(ni == sbi) {
			continue;
		}

		int nni = neigh.count[ni];

		// Subtract from mu the propensity of specie spi in subvolume
		// ni divided by nni (i.e. we sum a negative value).
		mu -= (rates.diffusion[GET_DR(spi, ni)]) / nni;

		// Add to sigma2 the propensity of specie spi in subvolume ni
		// divided by nni. No need to square since the coeff. is
		// always -1, just sum 1.
		sigma2 += (rates.diffusion[GET_DR(spi, ni)]) / nni;    
	}

	float m = max(EPSILON * x / g, 1.0f);
	float t1 = m / abs(mu);
	float t2 = (m * m) / (sigma2);

	return min(t1, t2);
}

__device__ float compute_tau_ncr(state state, reactions reactions, 
								 int * hors, bool crit_r[MAXREACTIONS], neigh neigh,
								 int sbi, rates rates)
{
	float min_tau = INFINITY;

	for (int spi = 0; spi < SPC; spi++) {

		// First of all we need to check if the specie spi is involved
		// in a critical event. If it is, skip it.
		bool skip_r = false;

		// check for critical reaction events
		for (int ri = 0; ri < RC; ri++) {
			// skip if reaction is critical and the specie is involved
			skip_r = skip_r || (crit_r[ri] && (reactions.r[GET_COEFF(spi, ri)] > 0));
		}

		// check for critical diffusion events
		bool skip_d = is_critical_diffusion(state, sbi, spi);

		if (skip_r && skip_d) { // ??? should be ||
			continue;
		}
		// spi is not involved in any critical event.

		float tau = compute_tau_sp(state, reactions, hors, crit_r, neigh, sbi, spi, rates);
		min_tau = min(min_tau, tau);
	}

	return min_tau;
}

__device__ float compute_tau_cr(state state, bool crit_r[MAXREACTIONS],
								int sbi, rates rates, curandStateMRG32k3a * s)
{
	float react_rates_sum_cr = 0.0;    // sum of the react rates of critical reactions
	float diff_rates_sum_cr = 0.0;     // sum of diffusion rates of critical diffusion events

	for (int ri = 0; ri < RC; ri++) {
		react_rates_sum_cr += (rates.reaction[GET_RR(ri, sbi)] * crit_r[ri]);
	}

	for (int spi = 0; spi < SPC; spi++) {
		diff_rates_sum_cr += (rates.diffusion[GET_DR(spi, sbi)] * is_critical_diffusion(state, sbi, spi));
	}

	if (react_rates_sum_cr == 0.0 && diff_rates_sum_cr == 0.0)
		return INFINITY;

	float rand = curand_uniform(&s[sbi]);
	return -logf(rand) / (react_rates_sum_cr + diff_rates_sum_cr);
}

__global__ void compute_taus(state state, reactions reactions, int * hors, neigh neigh,
							 rates rates, float * tau, float min_tau, char * leap, curandStateMRG32k3a * s)
{
	unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x;
	if (sbi >= SBC)
		return;

	// If on the previous step nobody changed our state from SSA_FF to
	// SSA, it means that no molecule entered here and we can
	// fast-forward time without recomputing anything.  Set tau =
	// old_tau - min_tau and return.
	//
	// Check:
	//     - If min_tau == 0.0 we are setting up the simulation, so pass.
	//     - If tau[sbi] is +Inf no event can happen here, so pass.
	//     - If min_tau == tau[sbi], it means that we have the same
	//       tau as the subvolume with the min_tau, but he was the one
	//       and we didn't get to act. We can't fast-forward (that
	//       would bring tau to zero), so just recompute a new tau.
	if (leap[sbi] == SSA_FF && !isinf(tau[sbi]) && min_tau > 0.0 && min_tau != tau[sbi]) {
		tau[sbi] -= min_tau;
		return;
	}

	// crit_r[ri] == TRUE if ri is critical in this subvolume
	bool crit_r[MAXREACTIONS]; 
	for (int ri = 0; ri < RC; ri++) {
		crit_r[ri] = is_critical_reaction(state, reactions, sbi, ri);
	}

	float tau_ncr = compute_tau_ncr(state, reactions, hors, crit_r, neigh, sbi, rates);
	float tau_cr = compute_tau_cr(state, crit_r, sbi, rates, s);

	// If tau_ncr is +Inf then every reaction is critical, and we
	// can't leap.  Also prevent leap if tau_ncr is too small.
	bool leap_here = true;
	if (isinf(tau_ncr) /*|| (tau_ncr < 10.0 / rates.matrix[GET_RATE(2, sbi)])*/) {
		// We start with fast-forward enabled. If someone diffuses to
		// us, they will need disable it by setting the state to SSA.
		leap[sbi] = SSA_FF;    
		leap_here = false;
	}

	if (tau_ncr < tau_cr) {
		// no critical event will happen, we'll leap with all the
		// non-critical events
		tau[sbi] = tau_ncr;
		if (leap_here)
			leap[sbi] = LEAP_NOCR;
	} else {
		// a single critical event will happen, all the non-critical
		// events will leap with tau = tau_cr
		tau[sbi] = tau_cr;
		if (leap_here)
			leap[sbi] = LEAP_CR;
	}

}