#ifndef CONSTANTS_CUH_
#define CONSTANTS_CUH_

#define MAXREACTIONS 64

extern __constant__ unsigned int SBC;
extern __constant__ int SPC;
extern __constant__ int RC;
extern __constant__ int NC;
extern __constant__ float EPSILON;

enum op: char {
	LEAP_CR,          // leap, then trigger a critical event
		LEAP_NOCR,    // leap
		SSA,          // plain SSA, tau is re-computed
		SSA_FF        // Fast-forwarded SSA, new_tau is old_tau - min_tau
		};


typedef struct s_rates {
	float * reaction;  // reaction  rates (rr)
	float * diffusion; // diffusion rates (dr)
	float * matrix;    // rates matrix {sum(rr), sum(dr), sum(rr) + sum(dr)}
	float * rc;        // reaction  constants
	float * dc;        // diffusion constants
} rates;

typedef struct s_reactions {
	int * r;   // reactants
	int * p;   // products
} reactions;


#endif /* CONSTANTS_CUH_ */
