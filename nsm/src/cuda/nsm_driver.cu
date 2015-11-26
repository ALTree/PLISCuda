#include "../../include/cuda/nsm_driver.cuh"

__constant__ unsigned int SBC;
__constant__ int SPC;
__constant__ int RC;
__constant__ int NC;
__constant__ float EPSILON;
__constant__ int * REACTANTS;

namespace NSMCuda {

// TODO: rename
// we keep nsm_step as a kernel, but this will be
// the main driver and "nsm" is not the right name.
void nsm(Topology t, State s, Reactions r, float * h_rrc, float * h_drc)
{
	unsigned int sbc = t.getN();
	int spc = s.getS();
	int rc = r.getR();

	int nc = 10;
	float epsilon = 0.05;

#if LOG
	std::cout << "\n   ***   Start simulation log   ***   \n\n";
#endif

	gpuErrchk(cudaMemcpyToSymbol(SBC, &sbc, sizeof(unsigned int)));
	gpuErrchk(cudaMemcpyToSymbol(SPC, &spc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(RC, &rc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(NC, &nc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(EPSILON, &epsilon, sizeof(float)));

#if LOG
	std::cout << "--- Allocating GPU memory... ";
#endif

	// ----- allocate and memcpy state array -----
	int * h_state = s.getState();

	int * d_state;
	gpuErrchk(cudaMalloc(&d_state, sbc * spc * sizeof(int)));
	gpuErrchk(cudaMemcpy(d_state, h_state, sbc * spc * sizeof(int), cudaMemcpyHostToDevice));

	// ----- allocate and memcpy reactants and products arrays -----
	int * h_reactants = r.getReactants();
	int * h_products = r.getProducts();

	int * d_reactants;
	int * d_products;
	gpuErrchk(cudaMalloc(&d_reactants, spc * rc * sizeof(int)));
	gpuErrchk(cudaMalloc(&d_products, spc * rc * sizeof(int)));
	gpuErrchk(cudaMemcpy(d_reactants, h_reactants, spc * rc * sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_products, h_products, spc * rc * sizeof(int), cudaMemcpyHostToDevice));

	// ----- allocate and memcpy topology array -----
	unsigned int * h_topology = t.getNeighboursArray();

	unsigned int * d_topology;
	gpuErrchk(cudaMalloc(&d_topology, 6 * sbc * sizeof(unsigned int)));
	gpuErrchk(cudaMemcpy(d_topology, h_topology, 6 * sbc * sizeof(unsigned int), cudaMemcpyHostToDevice));

	// ----- allocate rate matrix -----
	float * d_rate_matrix;
	gpuErrchk(cudaMalloc(&d_rate_matrix, 3 * sbc * sizeof(float)));

	// ----- allocate and memcpy rrc and drc -----
	float * d_rrc;
	float * d_drc;
	gpuErrchk(cudaMalloc(&d_rrc, rc * sizeof(float)));
	gpuErrchk(cudaMalloc(&d_drc, spc * sizeof(float)));
	gpuErrchk(cudaMemcpy(d_rrc, h_rrc, rc * sizeof(float), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_drc, h_drc, spc * sizeof(float), cudaMemcpyHostToDevice));

	// ----- allocate react_rates and diff_rates array
	float * d_react_rates_array;
	float * d_diff_rates_array;
	gpuErrchk(cudaMalloc(&d_react_rates_array, sbc * rc * sizeof(float)));
	gpuErrchk(cudaMalloc(&d_diff_rates_array, sbc * spc * sizeof(float)));

	// ----- allocate next_event thrust  vector
	thrust::device_vector<float> tau(sbc);

	// ----- allocate and initialize prng array
	curandStateMRG32k3a* d_prngstate;
	gpuErrchk(cudaMalloc(&d_prngstate, sbc * sizeof(curandStateMRG32k3a)));
	fill_prngstate_array<<<1, sbc>>>(d_prngstate);

	// ----- allocate leap array
	bool * d_leap;
	gpuErrchk(cudaMalloc(&d_leap, sbc * sizeof(bool)));

	// zero GPU memory, just to be sure
	// TODO: remove(?)
	gpuErrchk(cudaMemset(d_rate_matrix, 0, 3 * sbc * sizeof(float)));
	gpuErrchk(cudaMemset(d_react_rates_array, 0, sbc * rc * sizeof(float)));
	gpuErrchk(cudaMemset(d_diff_rates_array, 0, sbc * spc * sizeof(float)));
	gpuErrchk(cudaMemset(d_leap, 0, sbc * sizeof(bool)));

#if LOG
	std::cout << "done!\n";
	std::cout << "--- Initializing rate matrix... ";
#endif

	compute_rates<<<1, sbc>>>(d_state, d_reactants, d_topology, d_rate_matrix, d_rrc, d_drc, d_react_rates_array,
			d_diff_rates_array);

#if LOG
	std::cout << "done!\n";
#endif

	float * h_rate_matrix;

#if LOG
	h_rate_matrix = new float[3 * sbc];
	gpuErrchk(cudaMemcpy(h_rate_matrix, d_rate_matrix, 3 * sbc * sizeof(float), cudaMemcpyDeviceToHost));
	print_rate_matrix(h_rate_matrix, sbc);
#endif

#if LOG
	std::cout << "--- Fill initial next_event array... ";
#endif

	bool leap = true;

	if (!leap) {
		fill_tau_array<<<1, sbc>>>(thrust::raw_pointer_cast(tau.data()), d_rate_matrix);
	} else {
		fill_tau_array_leap<<<1, sbc>>>(d_state, d_reactants, d_products, d_topology, d_rate_matrix,
				d_react_rates_array, d_diff_rates_array, thrust::raw_pointer_cast(tau.data()), d_leap);
	}

#if LOG
	std::cout << "done!\n";
#endif

#if LOG
	print_tau(tau, sbc);
	std::cout << "--- Start simulation.\n\n";
#endif

	int steps = 4;

	for (int step = 1; step <= steps; step++) {

		int next = h_get_min_tau(tau);
		if (!leap) {
			nsm_step<<<1, sbc>>>(d_state, d_reactants, d_products, d_topology, d_rate_matrix, d_rrc, d_drc,
					d_react_rates_array, d_diff_rates_array, thrust::raw_pointer_cast(tau.data()), next, d_prngstate);
		} else {
			leap_step<<<1, sbc>>>(d_state, d_reactants, d_products, d_rate_matrix, d_topology, d_react_rates_array,
					d_diff_rates_array, d_rrc, d_drc, thrust::raw_pointer_cast(tau.data()), d_leap, d_prngstate);
		}

#if LOGSTEPS
		std::cout << "\n----- [step " << step << "] -----\n\n";

		// print system state
		gpuErrchk(cudaMemcpy(h_state, d_state, sbc * spc * sizeof(int), cudaMemcpyDeviceToHost));
		print_state(h_state, spc, sbc);

		// print rate matrix
		h_rate_matrix = new float[3 * sbc];
		gpuErrchk(cudaMemcpy(h_rate_matrix, d_rate_matrix, 3 * sbc * sizeof(float), cudaMemcpyDeviceToHost));
		print_rate_matrix(h_rate_matrix, sbc);

		// print tau array
		print_tau(tau, sbc);
#endif
	}
	gpuErrchk(cudaDeviceSynchronize());
#if LOG
	std::cout << "\n--- End simulation.\n\n";
#endif
}

// ----- print utils stuff -----

void print_state(int * h_state, int spc, int sbc)
{
	std::cout << "\n--- [system state] ---\n";
	for (int i = 0; i < sbc; i++) {
		std::cout << "sbv " << i << ": ";
		for (int j = 0; j < spc; j++)
			std::cout << h_state[j * sbc + i] << " ";
		std::cout << "\n";
	}
	std::cout << "----------------------\n";
}

void print_rate_matrix(float * h_rate_matrix, int sbc)
{
	std::cout << "\n--- [rate matrix] ---\n";
	for (int i = 0; i < sbc; i++) {
		std::cout << "sbv " << i << ": ";
		std::cout << h_rate_matrix[i] << " ";
		std::cout << h_rate_matrix[i + sbc] << " ";
		std::cout << h_rate_matrix[i + sbc * 2] << "\n";
	}
	std::cout << "---------------------\n\n";
}

void print_tau(thrust::device_vector<float> tau, int sbc)
{
	std::cout << "\n--- [tau array] ---\n";
	for (int i = 0; i < sbc; i++)
		std::cout << "sbv " << i << ": " << tau[i] << "\n";
	std::cout << "-------------------\n\n";
}

}
