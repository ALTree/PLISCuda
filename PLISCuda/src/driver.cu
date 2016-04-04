#include "../include/cuda/driver.cuh"

__constant__ unsigned int SBC;
__constant__ int SPC;
__constant__ int RC;
__constant__ int NC;
__constant__ float EPSILON;

namespace PLISCuda {

void run_simulation(Topology t, State s, Reactions r, float * h_rrc, float * h_drc, float endTime,
		int compartments_count, int * subv_constants, float log_freq)
{
	unsigned int sbc = t.getN();
	int spc = s.getS();
	int rc = r.getR();

	std::cout << "\n-- Starting Simulation -- \n\n";
	
	int threads = 512;
	int blocks = ceil(sbc / 512.0);
	std::cout << "  [using " << threads << " threads and " << blocks << " block(s)]\n\n";

	int nc = 10;    // threshold for critical/non-critical event
	float epsilon = 0.05;    // the epsilon parameter in the computation of the leap tau
							 // see [Cao06], formula 33

	std::cout << "  moving constants to Device...\n";

	// move constants to GPU
	gpuErrchk(cudaMemcpyToSymbol(SBC, &sbc, sizeof(unsigned int)));
	gpuErrchk(cudaMemcpyToSymbol(SPC, &spc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(RC, &rc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(NC, &nc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(EPSILON, &epsilon, sizeof(float)));

	std::cout << "  allocating memory on Device...\n";
	std::cout.precision(3);
	
	// ----- allocate and memcpy state arrays -----
	int * h_state = s.getState();

	int * d_state;
	gpuErrchk(cudaMalloc(&d_state, sbc * spc * sizeof(int)));
	gpuErrchk(cudaMemcpy(d_state, h_state, sbc * spc * sizeof(int), cudaMemcpyHostToDevice));

	int * d_state2;
	gpuErrchk(cudaMalloc(&d_state2, sbc * spc * sizeof(int)));

	std::cout << "    system state requires     " << (2*sbc*spc*sizeof(int)) / (1024.0 * 1024.0) << " MB\n";

	// ----- allocate and memcpy reactants and products arrays -----
	int * h_reactants = r.getReactants();
	int * h_products = r.getProducts();

	int * d_reactants;
	int * d_products;
	gpuErrchk(cudaMalloc(&d_reactants, spc * rc * sizeof(int)));
	gpuErrchk(cudaMalloc(&d_products, spc * rc * sizeof(int)));
	gpuErrchk(cudaMemcpy(d_reactants, h_reactants, spc * rc * sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_products, h_products, spc * rc * sizeof(int), cudaMemcpyHostToDevice));

	// ----- allocate and memcpy rrc and drc -----
	float * d_rrc;
	float * d_drc;
	gpuErrchk(cudaMalloc(&d_rrc, rc * sizeof(float) * compartments_count));
	gpuErrchk(cudaMalloc(&d_drc, spc * sizeof(float) * compartments_count));
	gpuErrchk(cudaMemcpy(d_rrc, h_rrc, rc * sizeof(float) * compartments_count, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_drc, h_drc, spc * sizeof(float) * compartments_count, cudaMemcpyHostToDevice));

	std::cout << "    reactions require         " <<
		(2*spc*rc*sizeof(int) + (rc+spc)*sizeof(float)) / (1024.0 * 1024.0) << " MB\n";
	
	// ----- allocate and memcpy topology array -----
	unsigned int * h_topology = t.getNeighboursArray();

	unsigned int * d_topology;
	gpuErrchk(cudaMalloc(&d_topology, 6 * sbc * sizeof(unsigned int)));
	gpuErrchk(cudaMemcpy(d_topology, h_topology, 6 * sbc * sizeof(unsigned int), cudaMemcpyHostToDevice));

	std::cout << "    system topology requires  " << (6*sbc*sizeof(unsigned int)) / (1024.0 * 1024.0) << " MB\n";
	
	// ----- allocate rate matrix -----
	float * d_rate_matrix;
	gpuErrchk(cudaMalloc(&d_rate_matrix, 3 * sbc * sizeof(float)));

	// ----- allocate and memcpy subv_constants array -----
	int * d_subv_consts;
	gpuErrchk(cudaMalloc(&d_subv_consts, sbc * sizeof(int)));
	gpuErrchk(cudaMemcpy(d_subv_consts, subv_constants, sbc * sizeof(int), cudaMemcpyHostToDevice));

	// ----- allocate react_rates and diff_rates array
	float * d_react_rates_array;
	float * d_diff_rates_array;
	gpuErrchk(cudaMalloc(&d_react_rates_array, sbc * rc * sizeof(float)));
	gpuErrchk(cudaMalloc(&d_diff_rates_array, sbc * spc * sizeof(float)));

	std::cout << "    rates arrays require      " << (sbc*(3+rc+spc)*sizeof(float)) / (1024.0 * 1024.0) << " MB\n";

	// ----- allocate tau thrust vector
	thrust::device_vector<float> tau(sbc);
	float * d_current_time;
	float h_current_time = 0.0;
	gpuErrchk(cudaMalloc(&d_current_time, sizeof(float)));

	// ----- allocate and initialize prng array
	curandStateMRG32k3a * d_prngstate;
	gpuErrchk(cudaMalloc(&d_prngstate, sbc * sizeof(curandStateMRG32k3a)));
	initialize_prngstate_array<<<blocks, threads>>>(d_prngstate);

	std::cout << "    prng array require        " << (sbc*sizeof(curandStateMRG32k3a)) / (1024.0 * 1024.0) << " MB\n";

	// ----- allocate leap and cr arrays
	char * d_leap;
	gpuErrchk(cudaMalloc(&d_leap, sbc * sizeof(char)));

	// zero GPU memory, just to be sure
	// TODO: remove(?) or check that we are zeroing everything
	gpuErrchk(cudaMemset(d_rate_matrix, 0, 3 * sbc * sizeof(float)));
	gpuErrchk(cudaMemset(d_react_rates_array, 0, sbc * rc * sizeof(float)));
	gpuErrchk(cudaMemset(d_diff_rates_array, 0, sbc * spc * sizeof(float)));
	gpuErrchk(cudaMemset(d_leap, 0, sbc * sizeof(char)));

	std::cout.precision(5);

	std::cout << "\n  computing initial rates...\n\n";
	
	compute_rates<<<blocks, threads>>>(d_state, d_reactants, d_topology, d_rate_matrix, d_rrc, d_drc, d_subv_consts,
			d_react_rates_array, d_diff_rates_array);

#ifdef DEBUG
	float * h_rate_matrix;
	h_rate_matrix = new float[3 * sbc];
	gpuErrchk(cudaMemcpy(h_rate_matrix, d_rate_matrix, 3 * sbc * sizeof(float), cudaMemcpyDeviceToHost));
	print_rate_matrix(h_rate_matrix, sbc);
#endif

	std::cout << "  computing initial taus...\n\n";

	fill_tau_array_leap<<<blocks, threads>>>(d_state, d_reactants, d_products, d_topology, d_rate_matrix,
			d_react_rates_array, d_diff_rates_array, thrust::raw_pointer_cast(tau.data()), 0.0, d_leap, d_prngstate);

#ifdef DEBUG
	print_tau(tau, sbc);
	print_leap_array(d_leap, sbc);
#endif

#ifndef LOG // we are in RELEASE mode
	std::ofstream out_file;
	float last_log_time = 0.0;
	if(log_freq > 0) { // we need to create the output file
		std::time_t t = std::time(0); // get timestamp to use in filename
		out_file.open("sim" + std::to_string(t) + ".data");
	}
#endif

	std::cout << "-- Begin Iterations -- \n\n";

	int step = 0;
	while(h_current_time < endTime) {

		step++;
		
#ifdef LOG
		std::cout << "  -- [step " << step << "] -- \n\n";
#endif

		// copy current state to d_state2
		gpuErrchk(cudaMemcpy(d_state2, d_state, spc * sbc * sizeof(int), cudaMemcpyDeviceToDevice));

		// get min tau from device
		int min_tau_sbi = h_get_min_tau(tau);
		if (isinf(tau[min_tau_sbi]) || tau[min_tau_sbi] < 0.0) {
			printf("\n -- WARNING: min(tau) = +Inf or < 0 - abort simulation --\n");
			break;
		}

		// update current time and forward it to the Device
		float min_tau = tau[min_tau_sbi];
		h_current_time += min_tau;
		gpuErrchk(cudaMemcpy(d_current_time, &h_current_time, sizeof(float), cudaMemcpyHostToDevice));

#ifdef LOG
		std::cout << "    tau = " << min_tau << "\n\n";
#endif

#ifndef LOG
		
#endif
		
		REPEAT:
		// first we leap, with tau = min_tau, in every subvolume that has leap enabled
		leap_step<<<blocks, threads>>>(d_state, d_reactants, d_products, d_rate_matrix, d_topology, d_react_rates_array,
				d_diff_rates_array, d_rrc, d_drc, tau[min_tau_sbi], d_current_time, d_leap, d_prngstate);

		// now we do a single ssa step, if min_tau was in a subvolume with leap not enabled
		ssa_step<<<blocks, threads>>>(d_state, d_reactants, d_products, d_topology, d_rate_matrix, d_react_rates_array,
				d_diff_rates_array, min_tau_sbi, d_current_time, d_leap, d_prngstate);

		// check if we need to revert this step
		thrust::device_vector<bool> revert(sbc);
		check_state<<<blocks, threads>>>(d_state, thrust::raw_pointer_cast(revert.data()));
		bool revert_state = !thrust::none_of(revert.begin(), revert.end(), thrust::identity<bool>());
		if(revert_state) {

#ifdef LOG
			std::cout << "\n -- WARNING: need to revert state --\n";
			std::cout << "\t old tau = " << min_tau << ", time was = " << h_current_time << "\n";
			std::cout << "\t new tau = " << min_tau / 2.0 << ", ";
#endif

			// restore state from the copy
			gpuErrchk(cudaMemcpy(d_state, d_state2, spc * sbc * sizeof(int), cudaMemcpyDeviceToDevice));

			// halven tau and update current time 
			h_current_time = h_current_time - min_tau + min_tau / 2.0;
			min_tau = min_tau / 2.0;

#ifdef LOG
			std::cout << "time is  = " << h_current_time << "\n";
#endif

			// send the new tau to the Device
			gpuErrchk(cudaMemcpy(d_current_time, &h_current_time, sizeof(float), cudaMemcpyHostToDevice));

			goto REPEAT;
		}

		// update rates
		// TODO: the computed values are not used if the subvolume
		// is tagged as SSA_FF, so we should avoid doing the
		// computation
		compute_rates<<<blocks, threads>>>(d_state, d_reactants, d_topology, d_rate_matrix, d_rrc, d_drc, d_subv_consts,
				d_react_rates_array, d_diff_rates_array);

		// update tau array
		fill_tau_array_leap<<<blocks, threads>>>(d_state, d_reactants, d_products, d_topology, d_rate_matrix,
				d_react_rates_array, d_diff_rates_array, thrust::raw_pointer_cast(tau.data()), min_tau, d_leap,
				d_prngstate);

#ifdef LOG
		gpuErrchk(cudaMemcpy(h_state, d_state, sbc * spc * sizeof(int), cudaMemcpyDeviceToHost));
		print_state(h_state, spc, sbc, h_current_time);
#endif
		
#ifdef DEBUG
		h_rate_matrix = new float[3 * sbc];
		gpuErrchk(cudaMemcpy(h_rate_matrix, d_rate_matrix, 3 * sbc * sizeof(float), cudaMemcpyDeviceToHost));
		print_rate_matrix(h_rate_matrix, sbc);

		print_tau(tau, sbc);
		print_leap_array(d_leap, sbc);
#endif

#ifndef LOG // log to file only when in RELEASE mode
		if(log_freq > 0 && (h_current_time - last_log_time) >= log_freq) {
			last_log_time = h_current_time;
			gpuErrchk(cudaMemcpy(h_state, d_state, sbc * spc * sizeof(int), cudaMemcpyDeviceToHost));
			std::string s = "";
			s += "-- state at " + std::to_string(h_current_time) + "\n";
			for (int i = 0; i < sbc; i++) {
				s += "sbv " + std::to_string(i) + ": ";
				for (int j = 0; j < spc; j++)
					s +=  std::to_string(h_state[j * sbc + i]) + " ";
				s += "\n";
			}
			s += "\n";
			out_file << s;
		}
#endif
		
	}

	out_file.close();
	
	gpuErrchk(cudaDeviceSynchronize());

	std::cout << "-- Simulation Complete -- \n";

	// print final state no matter which mode we are using
	gpuErrchk(cudaMemcpy(h_state, d_state, sbc * spc * sizeof(int), cudaMemcpyDeviceToHost));
	print_state(h_state, spc, sbc, h_current_time);

	std::cout << "  final simulation time: " << h_current_time << "\n";

}

int h_get_min_tau(thrust::device_vector<float> &tau)
	{
	thrust::device_vector<float>::iterator iter = thrust::min_element(tau.begin(), tau.end());
	return iter - tau.begin();
}

// ----- print utils stuff -----

void print_state(int * h_state, int spc, int sbc, float current_time)
{
	std::cout << "\n--- [state at t = " << current_time << "] ---\n";
	for (int i = 0; i < sbc; i++) {
		std::cout << "sbv " << i << ": ";
		for (int j = 0; j < spc; j++)
			std::cout << h_state[j * sbc + i] << " ";
		std::cout << "\n";
	}
	std::cout << "----------------------\n\n";
}

void print_rate_matrix(float * h_rate_matrix, int sbc)
{
	std::cout << "--- [rate matrix] ---\n";
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
	std::cout << "--- [tau array] ---\n";
	for (int i = 0; i < sbc; i++)
		std::cout << "sbv " << i << ": " << tau[i] << "\n";	
	std::cout << "-------------------\n\n";
}

void print_leap_array(char * d_leap, int sbc)
{
	std::cout << "--- [leap array] ---\n";
	char * h_leap = new char[sbc];
	std::vector<std::string> to_print = {"LEAP_CR", "LEAP_NOCR", "SSA", "SSA_FF"};
	gpuErrchk(cudaMemcpy(h_leap, d_leap, sbc * sizeof(char), cudaMemcpyDeviceToHost));
	for (int i = 0; i < sbc; i++) {
		std::cout << "sbi " << i << "] " << "leap: " <<
			to_print[h_leap[i] - LEAP_CR] << "\n";		
	}
	std::cout << "-------------------\n\n";
}

}
