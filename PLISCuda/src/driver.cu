#include "../include/cuda/driver.cuh"
#include <cuda_profiler_api.h>

__constant__ unsigned int SBC;
__constant__ int SPC;
__constant__ int RC;
__constant__ int NC;
__constant__ float EPSILON;

namespace PLISCuda {

	void run_simulation(Topology t, State s, Reactions r, float * h_rrc, float * h_drc, float endTime,
						int compartments_count, int * subv_constants, float log_freq)
	{
		// the number of subvolumes
		unsigned int sbc = t.getN();
		// sbc rounded up to the nearest multiple of 4 to shut up
		// initcheck about uninitialized memory accesses.
		unsigned int sbc4 = (sbc % 4 ? sbc+(4-sbc%4) : sbc);

		int spc = s.getS();   // number of species
		int rc = r.getR();    // number of reactions

		std::cout << "\n-- Starting Simulation -- \n\n";

		// We use 512 threads for block, except when we have just 1
		// block; in that case the number of threads will be equal to
		// the number of subvolumes rounded up to a multiple of 32.
		int threads = 512;
		int blocks = ceil(sbc / 512.0);
		if (blocks == 1)
			threads = ((sbc / 32) + 1)*32;

		std::cout << "  [using " << threads << " threads and " << blocks << " block(s)]\n\n";

		int nc = 10;             // threshold for critical/non-critical event
		float epsilon = 0.05;    // the epsilon parameter in the computation of the leap tau
		                         // see [Cao06], formula 33

		std::cout << "  moving constants to Device...\n";

		// move constants to the GPU
		gpuErrchk(cudaMemcpyToSymbol(SBC, &sbc, sizeof(unsigned int)));
		gpuErrchk(cudaMemcpyToSymbol(SPC, &spc, sizeof(int)));
		gpuErrchk(cudaMemcpyToSymbol(RC, &rc, sizeof(int)));
		gpuErrchk(cudaMemcpyToSymbol(NC, &nc, sizeof(int)));
		gpuErrchk(cudaMemcpyToSymbol(EPSILON, &epsilon, sizeof(float)));

		std::cout << "  allocating memory on Device...\n";
		std::cout.precision(3);


	
		// -------------------- state --------------------
		int * h_state = s.getState();
		int * d_state_curr;
		int * d_state_next;
		gpuErrchk(cudaMalloc(&d_state_curr, sbc4 * spc * sizeof(int)));
		gpuErrchk(cudaMalloc(&d_state_next, sbc4 * spc * sizeof(int)));

		// even if we allocated sbc4-long arrays, we must memcpy
		// *exactly* sbc values, or we'll hit uninitialized memory on
		// the CPU!
		gpuErrchk(cudaMemcpy(d_state_curr, h_state, sbc * spc * sizeof(int), cudaMemcpyHostToDevice));

		// state.next needs to be initialized to the old state
		gpuErrchk(cudaMemcpy(d_state_next, d_state_curr, sbc * spc * sizeof(int), cudaMemcpyDeviceToDevice));
		
		state state = {
			d_state_curr,
			d_state_next
		};

		std::cout << "    system state requires     " 
				  << (2*sbc*spc*sizeof(int))  // state (current and next)
			/ (1024.0 * 1024.0) 
				  << " MB\n";
		// -----------------------------------------------



		// -------------------- reactions --------------------
		// reactants and products
		int * h_reactants = r.getReactants();
		int * h_products = r.getProducts();
		int * d_reactants;
		int * d_products;
		gpuErrchk(cudaMalloc(&d_reactants, spc * rc * sizeof(int)));
		gpuErrchk(cudaMalloc(&d_products, spc * rc * sizeof(int)));
		gpuErrchk(cudaMemcpy(d_reactants, h_reactants, spc * rc * sizeof(int), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_products, h_products, spc * rc * sizeof(int), cudaMemcpyHostToDevice));

		reactions reactions = {
			d_reactants,
			d_products,
		};

		// reactions and diffusion constants
		float * d_rrc;
		float * d_drc;
		gpuErrchk(cudaMalloc(&d_rrc, rc * sizeof(float) * compartments_count));
		gpuErrchk(cudaMalloc(&d_drc, spc * sizeof(float) * compartments_count));
		gpuErrchk(cudaMemcpy(d_rrc, h_rrc, rc * sizeof(float) * compartments_count, cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_drc, h_drc, spc * sizeof(float) * compartments_count, cudaMemcpyHostToDevice));

		// hors cache
		int * d_hors;
		gpuErrchk(cudaMalloc(&d_hors, spc * sizeof(int)));

		// hors needs to be initialized immediately 
		init_hors<<<1, 1>>>(d_hors, reactions, spc);

		std::cout << "    reactions require         "
				  << (2*spc*rc*sizeof(int) +    // reactants and products
					  (rc+spc)*sizeof(float) +  // reactions and diffusion rates
					  (spc*sizeof(int)))        // hors
			/ (1024.0 * 1024.0)
				  << " MB\n";
		// -----------------------------------------------
	


		// -------------------- rates --------------------
		float * d_rate_matrix;
		float * d_react_rates_array;
		float * d_diff_rates_array;
		gpuErrchk(cudaMalloc(&d_rate_matrix, 3 * sbc * sizeof(float)));
		gpuErrchk(cudaMalloc(&d_react_rates_array, sbc * rc * sizeof(float)));
		gpuErrchk(cudaMalloc(&d_diff_rates_array, sbc * spc * sizeof(float)));
		
		rates rates = {
			d_react_rates_array, 
			d_diff_rates_array, 
			d_rate_matrix,
			d_rrc, d_drc
		};

		int * d_subv_consts;
		gpuErrchk(cudaMalloc(&d_subv_consts, sbc4 * sizeof(int)));
		gpuErrchk(cudaMemcpy(d_subv_consts, subv_constants, sbc * sizeof(int), cudaMemcpyHostToDevice));

		std::cout << "    rates arrays require      " 
				  << (sbc*(3+rc+spc)*sizeof(float)) 
			/ (1024.0 * 1024.0) 
				  << " MB\n";
		// -----------------------------------------------

		

		// -------------------- topology --------------------
		unsigned int * h_topology = t.getNeighboursArray();
		unsigned int * d_topology;
		int * d_neighcount;
		gpuErrchk(cudaMalloc(&d_topology, 6 * sbc * sizeof(unsigned int)));
		gpuErrchk(cudaMalloc(&d_neighcount, sbc * sizeof(int)));
		gpuErrchk(cudaMemcpy(d_topology, h_topology, 6 * sbc * sizeof(unsigned int), cudaMemcpyHostToDevice));

		neigh neigh = {
			d_topology,
			d_neighcount
		};

		// neigh.count needs to be initialized immediately
		init_ncount<<<blocks, threads>>>(neigh);

		std::cout << "    system topology requires  " << (sbc*sizeof(int) + 
														  6*sbc*sizeof(unsigned int)) / (1024.0 * 1024.0) 
				  << " MB\n";
		// --------------------------------------------------



		// -------------------- tau --------------------
		thrust::device_vector<float> tau(sbc4, INFINITY);
		float * d_current_time;
		float h_current_time = 0.0;
		gpuErrchk(cudaMalloc(&d_current_time, sizeof(float)));
		// ---------------------------------------------



		// -------------------- prng --------------------
		curandStateMRG32k3a * d_prngstate;
		gpuErrchk(cudaMalloc(&d_prngstate, sbc * sizeof(curandStateMRG32k3a)));
		init_prng<<<blocks, threads>>>(d_prngstate);

		std::cout << "    prng array require        " 
				  << (sbc*sizeof(curandStateMRG32k3a)) 
			/ (1024.0 * 1024.0) 
				  << " MB\n";
		// ----------------------------------------------



		// -------------------- leap --------------------
		char * d_leap;
		gpuErrchk(cudaMalloc(&d_leap, sbc4 * sizeof(char)));
		// -------------------- leap --------------------



		// -------------------- revert --------------------
		thrust::device_vector<int> revert(sbc4, 0);
		// ------------------------------------------------



		// zero GPU memory, just to be sure
		// TODO: remove(?) or check that we are zeroing everything
		gpuErrchk(cudaMemset(d_rate_matrix, 0, 3 * sbc * sizeof(float)));
		gpuErrchk(cudaMemset(d_react_rates_array, 0, sbc * rc * sizeof(float)));
		gpuErrchk(cudaMemset(d_diff_rates_array, 0, sbc * spc * sizeof(float)));
		gpuErrchk(cudaMemset(d_leap, SSA_FF, sbc4 * sizeof(char)));



		std::cout.precision(5);

		std::cout << "\n  computing initial rates...\n\n";

		compute_rates<<<blocks, threads>>>(state, reactions, neigh, rates, d_subv_consts);

#ifdef DEBUG
		float * h_rate_matrix;
		h_rate_matrix = new float[3 * sbc];
		gpuErrchk(cudaMemcpy(h_rate_matrix, d_rate_matrix, 3 * sbc * sizeof(float), cudaMemcpyDeviceToHost));
		print_rate_matrix(h_rate_matrix, sbc);
#endif

		std::cout << "  computing initial taus...\n\n";

		compute_taus<<<blocks, threads>>>(state, reactions, d_hors, neigh, rates, 
										  thrust::raw_pointer_cast(tau.data()), 0.0, 
										  d_leap, d_prngstate);

#ifdef DEBUG
		print_tau(tau, sbc);
		print_leap_array(d_leap, sbc);
#endif

#ifndef LOG // log to file only when in RELEASE mode 
		float last_log_time = 0.0;
		std::time_t tstamp = std::time(0); // get timestamp to use in filename
		std::clock_t last_progress_report = std::clock();
#endif

		std::cout << "-- Begin Iterations -- \n\n";

		std::clock_t sim_start = std::clock();


#ifdef PROFILE
		cudaProfilerStart();
#endif

		// ------------------------------ simulation loop start ------------------------------

		int step = 0;
		while(h_current_time < endTime) {
			step++;
		
#ifdef LOG
			std::cout << "  -- [step " << step << "] -- \n\n";
#endif

			// get min tau from device
			thrust::device_vector<float>::iterator iter = thrust::min_element(tau.begin(), tau.end());
			int min_tau_sbi = iter - tau.begin();
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

		REPEAT:
			// first we leap, with tau = min_tau, in every subvolume that has leap enabled
			leap_step<<<blocks, threads>>>(state, reactions, neigh, rates,
										   tau[min_tau_sbi], d_current_time, d_leap, d_prngstate);

			// now we do a single ssa step, if min_tau was in a subvolume with leap not enabled
			ssa_step<<<blocks, threads>>>(state, reactions, neigh, rates,
										  min_tau_sbi, d_current_time, d_leap, d_prngstate);

			// check if we need to revert this step
			check_state<<<blocks, threads>>>(state, thrust::raw_pointer_cast(revert.data()));
			bool revert_state = !thrust::none_of(revert.begin(), revert.end(), thrust::identity<int>());
			if(revert_state) {
#ifdef LOG
				std::cout << "\n -- WARNING: need to revert state --\n";
				std::cout << "\t old tau = " << min_tau << ", time was = " << h_current_time << "\n";
				std::cout << "\t new tau = " << min_tau / 2.0 << ", ";
#endif
				
				// reset state.next from curr, halven tau, and retry
				gpuErrchk(cudaMemcpy(state.next, state.curr, sbc * spc * sizeof(int), cudaMemcpyDeviceToDevice));
				h_current_time = h_current_time - min_tau + min_tau / 2.0;
				min_tau = min_tau / 2.0;

#ifdef LOG
				std::cout << "time is  = " << h_current_time << "\n";
#endif

				// send the new tau to the Device
				gpuErrchk(cudaMemcpy(d_current_time, &h_current_time, sizeof(float), cudaMemcpyHostToDevice));

				goto REPEAT;
			} // end if(revert_state)

			// If we're good (i.e. state.next has passed the
			// check_state call), swap state.next with state.curr, and
			// then copy state.curr to next
			int * temp = state.curr;
			state.curr = state.next;
			state.next = temp;
			gpuErrchk(cudaMemcpy(state.next, state.curr, sbc * spc * sizeof(int), cudaMemcpyDeviceToDevice));

			// update rates
			// TODO: the computed values are not used if the subvolume
			// is tagged as SSA_FF, so we should avoid doing the
			// computation
			compute_rates<<<blocks, threads>>>(state, reactions, neigh, rates, d_subv_consts);

			// update tau array
			compute_taus<<<blocks, threads>>>(state, reactions, d_hors, neigh, rates, 
											  thrust::raw_pointer_cast(tau.data()), min_tau, 
											  d_leap, d_prngstate);

#ifdef LOG
			gpuErrchk(cudaMemcpy(h_state, state.curr, sbc * spc * sizeof(int), cudaMemcpyDeviceToHost));
			print_state(h_state, spc, sbc, h_current_time);
#endif
		
#ifdef DEBUG
			h_rate_matrix = new float[3 * sbc];
			gpuErrchk(cudaMemcpy(h_rate_matrix, d_rate_matrix, 3 * sbc * sizeof(float), cudaMemcpyDeviceToHost));
			print_rate_matrix(h_rate_matrix, sbc);

			print_tau(tau, sbc);
			print_leap_array(d_leap, sbc);
#endif

#ifndef LOG			
			std::clock_t curr_time = std::clock();
			if (curr_time - last_progress_report > 30*CLOCKS_PER_SEC) {
				std::cout << "  [";
				print_eltime((curr_time - sim_start)/CLOCKS_PER_SEC);
				std::cout << "]  simulation time = ";
				printf("%6.3fs\n", h_current_time);
				last_progress_report = curr_time;
			} 
#endif


#ifndef LOG // log to file only when in RELEASE mode
			if(log_freq > 0 && (h_current_time - last_log_time) >= log_freq) {
				// get system state from the GPU
				last_log_time = h_current_time;
				gpuErrchk(cudaMemcpy(h_state, state.curr, sbc * spc * sizeof(int), cudaMemcpyDeviceToHost));
				
				// create a new file, named sim<timestamp>_<snapshot_time>.dat,
				// and write the system state in it
				std::ofstream log_file;
				log_file.open("sim" + std::to_string(tstamp) + "_" + std::to_string(h_current_time) +".dat");
				log_file << print_state_snapshot(h_state, spc, sbc, h_current_time);
				log_file.close();
			}
#endif
		}

		// ------------------------------ simulation loop end ------------------------------

#ifdef PROFILE
		cudaProfilerStop();
#endif

		gpuErrchk(cudaDeviceSynchronize());

		std::clock_t sim_end = std::clock();

		std::cout << "\n-- Simulation Complete -- \n\n";

#ifndef LOG // when logging to file, remember to print the final state
		gpuErrchk(cudaMemcpy(h_state, state.curr, sbc * spc * sizeof(int), cudaMemcpyDeviceToHost));

		std::ofstream log_file;
		log_file.open("sim" + std::to_string(tstamp) + "_" + std::to_string(h_current_time) +".dat");
		log_file << print_state_snapshot(h_state, spc, sbc, h_current_time);
		log_file.close();
#endif

		// print some final info to stdout
		std::cout << "  final simulation time:  " << h_current_time << "s\n";
		float eltime = float(sim_end - sim_start) / CLOCKS_PER_SEC;
		std::cout << "  elapsed time:           ";
		print_eltime(eltime);
		std::cout << "\n  simulation steps:       " << step << "\n";
		std::cout << "  steps/second:           " << step/eltime << "\n";

		std::cout << "\n";
	} 

	
	// -------------------- utils functions --------------------
  
	void print_eltime(float secs)
	{
		int hrs = secs / 3600;
		int mts = (secs - hrs*3600)/60;
		float scs = secs - (hrs*3600 + mts*60);
		std::cout << hrs << "h ";
		std::cout << mts << "m ";
		printf("%5.2fs", scs);
	}
	
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

	std::string print_state_snapshot(int * h_state, int spc, int sbc, float current_time)
	{
		std::string s = "";
		s += "subvolumes=" + std::to_string(sbc) + "\n";
		s += "species=" + std::to_string(spc) + "\n";
		s += "time=" + std::to_string(current_time) + "\n\n";
		for (int i = 0; i < sbc; i++) {
			for (int j = 0; j < spc; j++)
				s +=  std::to_string(h_state[j * sbc + i]) + " ";
			s += "\n";
		}
		return s;
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
 