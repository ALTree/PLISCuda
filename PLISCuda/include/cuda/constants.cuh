#ifndef CONSTANTS_CUH_
#define CONSTANTS_CUH_

extern __constant__ unsigned int SBC;
extern __constant__ int SPC;
extern __constant__ int RC;
extern __constant__ int NC;
extern __constant__ float EPSILON;

enum op: char {
	LEAP_CR,      // leap, then trigger a critical event
	LEAP_NOCR,    // leap
	SSA,          // plain SSA, tau is re-computed
	SSA_FF        // Fast-forwarded SSA, new_tau is old_tau - min_tau
};

#endif /* CONSTANTS_CUH_ */
