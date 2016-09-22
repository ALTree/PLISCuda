#include "../include/cuda/leap.cuh"

__global__ void leap_step(state state, reactions reactions, neigh neigh,
						  rates rates, float min_tau,
						  float * current_time, char * leap, curandStateMRG32k3a * prngstate)
{
	INDCHECK()

	if (leap[sbi] == SSA || leap[sbi] == SSA_FF)
		return;

	int neigh_count = neigh.count[sbi];

	// precompute is_critical_reaction
	bool crit_r[MAXREACTIONS];

	for (int ri = 0; ri < RC; ri++)
		crit_r[ri] = is_critical_reaction(state, reactions, sbi, ri);

	// fire all the non-critical reaction events
	for (int ri = 0; ri < RC; ri++) {

		if (crit_r[ri])
			continue;

		unsigned int k;    // how many times it fires
		k = curand_poisson(&prngstate[sbi], min_tau * rates.reaction[GET_RR(ri, sbi)]);

#ifdef LOG
		if(k > 0)
			printf("(%f) [subv %d] fire reaction %d for %d times\n", *current_time, sbi, ri, k);
#endif
		
		// update state
		for (int spi = 0; spi < SPC; spi++) {
			int new_state = k * (reactions.p[GET_COEFF(spi, ri)] - reactions.r[GET_COEFF(spi, ri)]);
			atomicAdd(&state.next[GET_SPI(spi, sbi)], new_state);
		}

	}

	__syncthreads();

	// fire non-critical outgoing diffusion events
	for (int spi = 0; spi < SPC; spi++) {
		if (is_critical_diffusion(state, sbi, spi))
			continue;

		// diffuse to each neighbour
		for (unsigned int ngb = 0; ngb < neigh_count; ngb++) {
			unsigned int k = curand_poisson(&prngstate[sbi], min_tau * rates.diffusion[GET_DR(spi, sbi)]);
			atomicAdd(&state.next[GET_SPI(spi, neigh.index[sbi*6 + ngb])], k);
			atomicSub(&state.next[GET_SPI(spi, sbi)], k);
			if (leap[neigh.index[sbi * 6 + ngb]] == SSA_FF) {
				leap[neigh.index[sbi * 6 + ngb]] = SSA;    // set the OP of the receiver to SSA
			}

#ifdef LOG
			if (k > 0)
				printf("(%f) [subv %d] diffuse %d molecules of specie %d to subv %d\n",
					   *current_time, sbi, k, spi, neigh.index[sbi * 6 + ngb]);
#endif
		}

	}

	if (leap[sbi] != LEAP_CR)
		return;

	// TODO: this comment is outdated
	// Problem: in the following we use the old react_rates_array
	// (we'll update it later, in another kernel call), but it can
	// happen that a reaction has rate > 0 even if we can't fire it
	// (for example because we had 2 molecules at the start of this
	// step, but the leap diffusion phase we just executed removed one
	// of them).
	// 
	// To avoid negative population, we choose the random reaction as
	// usual, but we only fire it if we have enough molecules in the
	// subvolume.  Note that it's possibile that something has entered
	// from a neighbour during the leap phase, so checking the state
	// is a good pragmatic way to ensure that we'll fire everytime is
	// possibile.

	__syncthreads();

	// fire a single critical event
	float rand = curand_uniform(&prngstate[sbi]);

	// first we have to choose if we'll fire a reaction or we'll
	// diffuse a molecule

	// sum the reaction rates of critical reactions
	float rr_sum = 0.0;
	for (int ri = 0; ri < RC; ri++)
		rr_sum += rates.reaction[GET_RR(ri, sbi)] * crit_r[ri];

	// sum the diffusion rates of critical diffusion events
	float dr_sum = 0.0;
	for (int spi = 0; spi < SPC; spi++)
		dr_sum += rates.diffusion[GET_DR(spi, sbi)] * is_critical_diffusion(state, sbi, spi);

	// if the sum is zero we can't fire or diffuse anything
	if (rr_sum + dr_sum == 0.0) {
		return;
	}

	if (rand < rr_sum / (rr_sum + dr_sum)) {    // reaction

		float scaled_sum = rr_sum * rand;
		float partial_sum = 0;

		int ric = 0;
		while (partial_sum <= scaled_sum) {
			partial_sum += rates.reaction[GET_RR(ric, sbi)] * crit_r[ric];
			ric++;
		}
		// We'll fire the ric-nth critical reactions.

		int ri;
		for (ri = 0; ri < RC && ric > 0; ri++) {
			if (crit_r[ri])
				ric--;
		}
		ri = ri - 1;

		// Check if the current state lets us fire reaction ri (see
		// comment above).
		bool fire = true;
		for (int spi = 0; spi < SPC; spi++) {
			fire = fire && (state.curr[GET_SPI(spi, sbi)] >= reactions.r[GET_COEFF(spi, ri)]);
		}

		if (!fire)
			return;

		for (int spi = 0; spi < SPC; spi++) {
			int delta = reactions.p[GET_COEFF(spi, ri)] - reactions.r[GET_COEFF(spi, ri)];
			atomicAdd(&state.next[GET_SPI(spi, sbi)], delta);
		}

#ifdef LOG
		printf("(%f) [subv %d] fire reaction %d (critical)\n", *current_time, sbi, ri);
#endif

	} else {    // diffusion

		float scaled_sum = dr_sum * rand;    // no need to scale down by neigh_count (it's the raw sum)
		float partial_sum = 0;

		int spic = 0;
		while (partial_sum <= scaled_sum) {
			partial_sum += rates.diffusion[GET_DR(spic, sbi)] * is_critical_diffusion(state, sbi, spic);
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

		// Check if the current state lets us diffuse specie spi (see
		// comment above).
		bool fire = state.curr[GET_SPI(spi, sbi)] > 0;

		if (!fire)
			return;

		// Choose a random destination.
		// We should re-use the rand we already have, but it doesn't
		// really matter, since we are using rejection sampling to
		// choose the random destination.
		int rdi;
		do {
			rdi = (int) (curand_uniform(&prngstate[sbi]) * 6);
		} while (rdi > 5);

		// get index of neighbour #rdi (overwrite rdi, whatever)
		rdi = neigh.index[sbi * 6 + rdi];

		// If rdi == sbi (i.e. diffuse to myself) don't do anything
		if (rdi != sbi) {
			atomicSub(&state.next[GET_SPI(spi, sbi)], 1);
			atomicAdd(&state.next[GET_SPI(spi, rdi)], 1);
			if (leap[neigh.index[rdi]] == SSA_FF) // TODO: atomic?
				leap[neigh.index[rdi]] = SSA;    // set the OP of the receiver to SSA
		}

#ifdef LOG
		printf("(%f) [subv %d] diffuse specie %d to %d (critical)\n", *current_time, sbi, spi, rdi);
#endif
	}

}

__global__ void check_state(state state, int * revert)
{
	INDCHECK()

	bool _revert = false;
	for (int spi = 0; spi < SPC; spi++)
		_revert = _revert || (state.next[GET_SPI(spi, sbi)] < 0);

	revert[sbi] = _revert ? 1 : 0;
}


