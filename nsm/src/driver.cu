#include "../../include/cuda/driver.cuh"

__constant__ unsigned int SBC;
__constant__ int SPC;
__constant__ int RC;
__constant__ int NC;
__constant__ float EPSILON;
__constant__ int * REACTANTS;
__constant__ bool LOG_EVENTS;

namespace NSMCuda {

void run_simulation(Topology t, State s, Reactions r, float * h_rrc, float * h_drc, int steps,
		int constants_files_count, int * subv_constants, struct ToLog to_log)
{
	unsigned int sbc = t.getN();
	int spc = s.getS();
	int rc = r.getR();

	int threads = 512;
	int blocks = ceil(sbc / 512.0);
	std::cout << "threads: " << threads << ", blocks: " << blocks << "\n";

	std::cout << "Will log data from subvolumes:\n\t";
	for (int i = 0; i < to_log.subv_len; i++)
		std::cout << to_log.subv[i] << " ";
	std::cout << "\n";

	std::cout << "and species\n\t";
	for (int i = 0; i < spc; i++)
		std::cout << (to_log.spc[i] ? std::to_string(i) : "") << " ";
	std::cout << "\n";

	std::cout << "with frequency\n\t";
	std::cout << to_log.freq << "\n\n";

	int nc = 10;    // threshold for critical/non-critical event
	float epsilon = 0.05;    // the epsilon parameter in the computation of the leap tau
							 // see [Cao06], formula 33

	bool log_events = false;

#if LOG
	std::cout << "\n   ***   Start simulation log   ***   \n\n";
#endif

	// move constants to GPU
	gpuErrchk(cudaMemcpyToSymbol(SBC, &sbc, sizeof(unsigned int)));
	gpuErrchk(cudaMemcpyToSymbol(SPC, &spc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(RC, &rc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(NC, &nc, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(EPSILON, &epsilon, sizeof(float)));
	gpuErrchk(cudaMemcpyToSymbol(LOG_EVENTS, &log_events, sizeof(bool)));

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
	gpuErrchk(cudaMalloc(&d_rrc, rc * sizeof(float) * constants_files_count));
	gpuErrchk(cudaMalloc(&d_drc, spc * sizeof(float) * constants_files_count));
	gpuErrchk(cudaMemcpy(d_rrc, h_rrc, rc * sizeof(float) * constants_files_count, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_drc, h_drc, spc * sizeof(float) * constants_files_count, cudaMemcpyHostToDevice));

	// ----- allocate and memcpy subv_constants array -----
	int * d_subv_consts;
	gpuErrchk(cudaMalloc(&d_subv_consts, sbc * sizeof(int)));
	gpuErrchk(cudaMemcpy(d_subv_consts, subv_constants, sbc * sizeof(int), cudaMemcpyHostToDevice));

	// ----- allocate react_rates and diff_rates array
	float * d_react_rates_array;
	float * d_diff_rates_array;
	gpuErrchk(cudaMalloc(&d_react_rates_array, sbc * rc * sizeof(float)));
	gpuErrchk(cudaMalloc(&d_diff_rates_array, sbc * spc * sizeof(float)));

	// ----- allocate tau thrust vector
	thrust::device_vector<float> tau(sbc);
	float * d_current_time;
	float h_current_time = 0.0;
	gpuErrchk(cudaMalloc(&d_current_time, sizeof(float)));

	// ----- allocate and initialize prng array
	curandStateMRG32k3a * d_prngstate;
	gpuErrchk(cudaMalloc(&d_prngstate, sbc * sizeof(curandStateMRG32k3a)));
	initialize_prngstate_array<<<blocks, threads>>>(d_prngstate);

	// ----- allocate leap and cr arrays
	char * d_leap;
	gpuErrchk(cudaMalloc(&d_leap, sbc * sizeof(char)));

	// ----- allocate and initialize log arrays
	unsigned int * d_sbv_to_log;
	bool * d_spc_to_log;
	gpuErrchk(cudaMalloc(&d_sbv_to_log, to_log.subv_len * sizeof(unsigned int)));
	gpuErrchk(cudaMalloc(&d_spc_to_log, spc * sizeof(bool)));
	gpuErrchk(cudaMemcpy(d_sbv_to_log, to_log.subv, to_log.subv_len * sizeof(unsigned int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_spc_to_log, to_log.spc, spc * sizeof(bool), cudaMemcpyHostToDevice));

	int MAX_LOG = 1000;
	int * d_log_data;
	gpuErrchk(cudaMalloc(&d_log_data, to_log.subv_len * to_log.spc_len * MAX_LOG * sizeof(int)));

	// zero GPU memory, just to be sure
	// TODO: remove(?) or check that we are zeroing everything
	gpuErrchk(cudaMemset(d_rate_matrix, 0, 3 * sbc * sizeof(float)));
	gpuErrchk(cudaMemset(d_react_rates_array, 0, sbc * rc * sizeof(float)));
	gpuErrchk(cudaMemset(d_diff_rates_array, 0, sbc * spc * sizeof(float)));
	gpuErrchk(cudaMemset(d_leap, 0, sbc * sizeof(char)));
	gpuErrchk(cudaMemset(d_log_data, 0, to_log.subv_len * to_log.spc_len * MAX_LOG * sizeof(int)));

	int * d_log_data_start = d_log_data;

#if LOG
	std::cout << "done!\n";
	std::cout << "--- Initializing rate matrix... ";
#endif

	compute_rates<<<blocks, threads>>>(d_state, d_reactants, d_topology, d_rate_matrix, d_rrc, d_drc, d_subv_consts,
			d_react_rates_array, d_diff_rates_array);

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

	fill_tau_array_leap<<<blocks, threads>>>(d_state, d_reactants, d_products, d_topology, d_rate_matrix,
			d_react_rates_array, d_diff_rates_array, thrust::raw_pointer_cast(tau.data()), 0.0, d_leap, d_prngstate);

#if LOG
	print_tau(tau, sbc);

	char * h_leap = new char[sbc];
	gpuErrchk(cudaMemcpy(h_leap, d_leap, sbc * sizeof(char), cudaMemcpyDeviceToHost));
	for (int i = 0; i < sbc; i++) {
		std::cout << "sbi " << i << "] " << "leap: " << h_leap[i] << "\n";
	}
#endif

#if LOG
	std::cout << "done!\n";
#endif

#if LOG
	print_tau(tau, sbc);
	std::cout << "--- Start simulation.\n\n";
#endif

	int log_counter = 0;
	float time_since_last_log = 0.0;
	for (int step = 1; step <= steps; step++) {

		if(h_current_time > 1000.0) {
			break;
		}

#if LOGSTEPS
		std::cout << "\n----- [step " << step << "] -----\n\n";
#endif

		if (false && !log_events && step % 100 == 0) {
			if (step > 100) {
				for (int i = 0; i < 15; i++)
					std::cout << "\b \b";
			}
			std::cout << std::setw(7) << step << "/" << steps;

		}

		// copy current state to d_state2
		gpuErrchk(cudaMemcpy(d_state2, d_state, spc * sbc * sizeof(int), cudaMemcpyDeviceToDevice));

		// get min tau from device
		int min_tau_sbi = h_get_min_tau(tau);
		if (isinf(tau[min_tau_sbi]) || tau[min_tau_sbi] < 0.0) {
			printf("\n\n--------------- WARNING: min(tau) = +Inf or < 0 - abort simulation ---------------\n\n");
			break;
		}

		// update current time and forward it to the Device
		float min_tau = tau[min_tau_sbi];
		h_current_time += min_tau;
		gpuErrchk(cudaMemcpy(d_current_time, &h_current_time, sizeof(float), cudaMemcpyHostToDevice));

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
		if (revert_state) {
#if LOGSTEPS
			std::cout << "\n--------------- REVERT STATE ---------------\n\n";
			std::cout << "----- old tau = " << min_tau << "time was = " << h_current_time << "\n";
			std::cout << "----- new tau = " << min_tau / 2.0 << " ";
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

		if (time_since_last_log >= to_log.freq && log_counter < MAX_LOG) {
			log_data<<<blocks, threads>>>(d_state, d_sbv_to_log, to_log.subv_len, d_spc_to_log, to_log.spc_len,
					d_log_data);
			time_since_last_log = 0.0;
			d_log_data = &d_log_data[to_log.spc_len * to_log.subv_len];
			log_counter++;
		} else {
			time_since_last_log += min_tau;
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

#if LOG
		std::cout << "\n";
#endif

#if LOGSTEPS

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
		gpuErrchk(cudaMemcpy(h_leap, d_leap, sbc * sizeof(bool), cudaMemcpyDeviceToHost));
		for (int i = 0; i < sbc; i++) {
			std::cout << "sbi " << i << "] " << "leap: " << h_leap[i] << "\n";
		}

#endif

	}

	gpuErrchk(cudaDeviceSynchronize());

#if LOG
	std::cout << "\n--- End simulation.\n\n";
#endif

	gpuErrchk(cudaMemcpy(h_state, d_state, sbc * spc * sizeof(int), cudaMemcpyDeviceToHost));
	std::cout << "\n";
	std::cout << "--- Final State ---\n";
	print_state(h_state, spc, sbc);
	std::cout << "Final simulation time: " << h_current_time << "\n";

	if (log_events) {
		std::cout << "\n\n--- Log Data ---\n";
		int * h_log_data = new int[log_counter * to_log.spc_len * to_log.subv_len];
		gpuErrchk(
				cudaMemcpy(h_log_data, d_log_data_start, log_counter * to_log.spc_len * to_log.subv_len * sizeof(int),
						cudaMemcpyDeviceToHost));

		for (int i = 0; i < 100; i++)
			std::cout << h_log_data[i] << " ";
		std::cout << "\n";

		for (int t = 0; t < log_counter; t++) {
			std::cout << "--- time ~ " << t * to_log.freq << "\n";
			for (int spi = 0; spi < to_log.spc_len; spi++) {
				std::cout << "\tlogged specie " << spi << "\n\t\t";
				for (int sbi = 0; sbi < to_log.subv_len; sbi++) {
					std::cout << h_log_data[spi * to_log.spc_len + sbi] << " ";
				}
				std::cout << "\n";

			}
			h_log_data = &h_log_data[to_log.spc_len * to_log.subv_len];
			std::cout << "\n";
		}
	}

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
