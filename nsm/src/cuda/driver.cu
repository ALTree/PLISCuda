#include "../../include/cuda/driver.cuh"

__constant__ unsigned int SBC;
__constant__ int SPC;
__constant__ int RC;
__constant__ int NC;
__constant__ float EPSILON;
__constant__ int * REACTANTS;
float TAU;

namespace NSMCuda {

void run_simulation(Topology t, State s, Reactions r, float * h_rrc, float * h_drc, int steps)
{
	unsigned int sbc = t.getN();
	int spc = s.getS();
	int rc = r.getR();

	int nc = 10;    // threshold for critical/non-critical event
	float epsilon = 0.05;    // the epsilon parameter in the computation of the leap tau
							 // see [Tao06], formula 33

#if LOG
	std::cout << "\n   ***   Start simulation log   ***   \n\n";
#endif

	// move constants to GPU
	gpuErrchk(cudaMemcpyToSymbol(SBC, &sbc, sizeof(unsigned int)));
	gpuErrchk(cudaMemcpyToSymbol(SPC, &spc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(RC, &rc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(NC, &nc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(EPSILON, &epsilon, sizeof(float)));

#if LOG
	std::cout << "--- Allocating GPU memory... ";
#endif

	// ----- allocate and memcpy state arrays -----
	int * h_state = s.getState();

	int * d_state;
	gpuErrchk(cudaMalloc(&d_state, sbc * spc * sizeof(int)));
	gpuErrchk(cudaMemcpy(d_state, h_state, sbc * spc * sizeof(int), cudaMemcpyHostToDevice));

	int * d_state2;
	gpuErrchk(cudaMalloc(&d_state2, sbc * spc * sizeof(int)));

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

	// ----- allocate tau thrust vector and current_time
	thrust::device_vector<float> tau(sbc);
	float * d_current_time;
	float h_current_time = 0.0;
	gpuErrchk(cudaMalloc(&d_current_time, sizeof(float)));

	// ----- allocate and initialize prng array
	curandStateMRG32k3a * d_prngstate;
	gpuErrchk(cudaMalloc(&d_prngstate, sbc * sizeof(curandStateMRG32k3a)));
	initialize_prngstate_array<<<1, sbc>>>(d_prngstate);

	// ----- allocate leap and cr arrays
	bool * d_leap;
	bool * d_cr;
	gpuErrchk(cudaMalloc(&d_leap, sbc * sizeof(bool)));
	gpuErrchk(cudaMalloc(&d_cr, sbc * sizeof(bool)));

	// zero GPU memory, just to be sure
	// TODO: remove(?) or check that we are zeroing everything
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

#if LOG
	float * h_rate_matrix;
	h_rate_matrix = new float[3 * sbc];
	gpuErrchk(cudaMemcpy(h_rate_matrix, d_rate_matrix, 3 * sbc * sizeof(float), cudaMemcpyDeviceToHost));
	print_rate_matrix(h_rate_matrix, sbc);
#endif

#if LOG
	std::cout << "--- Fill initial next_event array... ";
#endif

	fill_tau_array_leap<<<1, sbc>>>(d_state, d_reactants, d_products, d_topology, d_rate_matrix, d_react_rates_array,
			d_diff_rates_array, thrust::raw_pointer_cast(tau.data()), d_leap, d_cr, d_prngstate);

#if LOG
	print_tau(tau, sbc);

	bool * h_leap = new bool[sbc];
	bool * h_cr = new bool[sbc];
	gpuErrchk(cudaMemcpy(h_leap, d_leap, sbc * sizeof(bool), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(h_cr, d_cr, sbc * sizeof(bool), cudaMemcpyDeviceToHost));
	for (int i = 0; i < sbc; i++) {
		std::cout << "sbi " << i << "] " << "leap: " << (h_leap[i] ? "yes" : "no") << ", cr: "
				<< (h_cr[i] ? "yes" : "no") << "\n";
	}
#endif

#if LOG
	std::cout << "done!\n";
#endif

#if LOG
	print_tau(tau, sbc);
	std::cout << "--- Start simulation.\n\n";
#endif

	for (int step = 1; step <= steps; step++) {

		// copy current state to d_state2
		gpuErrchk(cudaMemcpy(d_state2, d_state, spc * sbc * sizeof(int), cudaMemcpyDeviceToDevice));

		// get min tau from device
		int min_tau_sbi = h_get_min_tau(tau);
		if (isinf(tau[min_tau_sbi])) {
			printf("\n\n--------------- WARNING: min(tau) = +Inf - abort simulation ---------------\n\n");
			break;
		}

		// update current time and forward it to the Device
		float min_tau = tau[min_tau_sbi];
		h_current_time += min_tau;
		gpuErrchk(cudaMemcpy(d_current_time, &h_current_time, sizeof(float), cudaMemcpyHostToDevice));

		REPEAT:
		// first we leap, with tau = min_tau, in every subvolume that has leap enabled
		leap_step<<<1, sbc>>>(d_state, d_reactants, d_products, d_rate_matrix, d_topology, d_react_rates_array,
				d_diff_rates_array, d_rrc, d_drc, tau[min_tau_sbi], d_current_time, d_leap, d_cr, d_prngstate);

		// now we do a single ssa step, if min_tau was in a subvolume with leap not enabled
		ssa_step<<<1, sbc>>>(d_state, d_reactants, d_products, d_topology, d_rate_matrix, d_react_rates_array,
				d_diff_rates_array, min_tau_sbi, d_current_time, d_leap, d_prngstate);

		// check if we need to revert this step
		thrust::device_vector<bool> revert(sbc);
		check_state<<<1, sbc>>>(d_state, thrust::raw_pointer_cast(revert.data()));
		bool revert_state = !thrust::none_of(revert.begin(), revert.end(), thrust::identity<bool>());
		if (revert_state) {
#if LOGSTEPS
			std::cout << "\n--------------- REVERT STATE ---------------\n\n";
			std::cout << "----- old tau = " << min_tau << "time was = " << h_current_time << "\n";
			std::cout << "----- new tau = " << min_tau/2.0 << " ";
#endif
			// restore state from the copy
			gpuErrchk(cudaMemcpy(d_state, d_state2, spc * sbc * sizeof(int), cudaMemcpyDeviceToDevice));

			h_current_time = h_current_time - min_tau + min_tau / 2.0;
			min_tau = min_tau / 2.0;
#if LOGSTEPS
			std::cout << "time is  = " << h_current_time << "\n";
			gpuErrchk(cudaMemcpy(d_current_time, &h_current_time, sizeof(float), cudaMemcpyHostToDevice));
#endif
			goto REPEAT;
		}

		// update rates
		compute_rates<<<1, sbc>>>(d_state, d_reactants, d_topology, d_rate_matrix, d_rrc, d_drc, d_react_rates_array,
				d_diff_rates_array);

		// update tau array
		fill_tau_array_leap<<<1, sbc>>>(d_state, d_reactants, d_products, d_topology, d_rate_matrix,
				d_react_rates_array, d_diff_rates_array, thrust::raw_pointer_cast(tau.data()), d_leap, d_cr,
				d_prngstate);

#if LOG
		std::cout << "\n";
#endif

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

		// print cr and leap arrays
		bool * h_leap = new bool[sbc];
		bool * h_cr = new bool[sbc];
		gpuErrchk(cudaMemcpy(h_leap, d_leap, sbc * sizeof(bool), cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(h_cr, d_cr, sbc * sizeof(bool), cudaMemcpyDeviceToHost));
		for (int i = 0; i < sbc; i++) {
			std::cout << "sbi " << i << "] " << "leap: " << (h_leap[i] ? "yes" : "no") << ", cr: "
			<< (h_cr[i] ? "yes" : "no") << "\n";
		}

#endif

	}

	gpuErrchk(cudaDeviceSynchronize());

#if LOG
	std::cout << "\n--- End simulation.\n\n";
#endif

	gpuErrchk(cudaMemcpy(h_state, d_state, sbc * spc * sizeof(int), cudaMemcpyDeviceToHost));
	print_state(h_state, spc, sbc);
}

int h_get_min_tau(thrust::device_vector<float> &tau)
{
	thrust::device_vector<float>::iterator iter = thrust::min_element(tau.begin(), tau.end());
	return iter - tau.begin();
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
