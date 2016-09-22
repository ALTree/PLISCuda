#ifndef CONSTANTS_CUH_
#define CONSTANTS_CUH_

#define INDCHECK() \
		unsigned int sbi = blockIdx.x * blockDim.x + threadIdx.x; \
		if (sbi >= SBC) return; 

// ---- model limitations ----
#define MAXREACTIONS 64


// ---- global device constants ----

extern __constant__ unsigned int SBC;
extern __constant__ int SPC;
extern __constant__ int RC;
extern __constant__ int NC;
extern __constant__ float EPSILON;


// ---- accessors macros ----

// Get specie_count index from state using specie index, subvolume
// index and subvolume count.
//
//     state[CUDA_GET_SPI(spi, sbi)];
#define GET_SPI(spi, sbi) ((spi) * (SBC) + (sbi))

// Get R or D or R+D index from rate matrix using subvolume index and
// 0, 1, 2 (resp. R, D, R+D).
//
//     rate_matrix[GET_RATE(i, sbi)]
#define GET_RATE(i, sbi) ((i) * (SBC) + (sbi))

// Get the stechiometric coefficient for specie spi in reaction ri.
//
//     reactants[GET_COEFF(spi, ri)]
#define GET_COEFF(spi, ri) ((spi) * (RC) + (ri))

// Get the react rate for reaction ri in subvolume sbi.
// 
//     react_rates_array[(GET_RR(ri, sbi)]
#define GET_RR(ri, sbi) ((ri) * (SBC) + (sbi))

// Get the diff rate for specie spi in subvolume sbi
// 
//     diff_rates_array[(GET_DR(spi, sbi)]
#define GET_DR(spi, sbi) ((spi) * (SBC) + (sbi))


// ---- structs for system data ----

enum op: char {
	LEAP_CR,          // leap, then trigger a critical event
		LEAP_NOCR,    // leap
		SSA,          // plain SSA, tau is re-computed
		SSA_FF        // Fast-forwarded SSA, new_tau is old_tau - min_tau
		};


typedef struct s_state {
	int * curr;   // current (previous step) system state
	int * next;   // next state (to be computed)
} state;

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

typedef struct s_neigh {
	unsigned int * index;   // neighboours indices 
	int * count;            // neighboours count
} neigh;


struct ToLog {
	unsigned int * subv;
	int subv_len;
	bool * spc;
	int spc_len;
	float freq;
};

#endif /* CONSTANTS_CUH_ */
