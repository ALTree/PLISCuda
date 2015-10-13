#include <cuda_runtime.h>

#include "../../include/cuda/nsm_driver.hpp"
#include "../../include/cuda/rates.cuh"
#include "../../include/cuda/cuda_utils.cuh"

namespace NSMCuda {

void nsm(Topology t, State s, Reactions r, float * reaction_rates_constants, float * diffusion_rates_constants)
{

	int reactions_count = r.getR();
	int subvolumes_count = t.getN();
	int species_count = s.getS();

	std::cout << "----- Allocating GPU memory ...";

	// ----- allocate state array -----
	int * h_state_array = s.getState();
	int state_len = s.getN() * s.getS();

	int * d_state_array;
	gpuErrchk(cudaMalloc(&d_state_array, state_len));

	// ----- allocate reactants and products arrays -----
	int * h_reactants_array = r.getReactants();
	int * h_products_array = r.getProducts();
	int reaction_len = r.getR() * r.getS();

	int * d_reactants_array, *d_products_array;
	gpuErrchk(cudaMalloc(&d_reactants_array, reaction_len));
	gpuErrchk(cudaMalloc(&d_products_array, reaction_len));

	// ----- allocate topology arrays -----
	int * h_topology_array = t.getNeighboursArray();
	int topology_len = 6 * t.getN();

	int * d_topology_array;
	gpuErrchk(cudaMalloc(&d_topology_array, topology_len));
	// ----- allocate rate matrix -----
	float * d_rate_matrix;
	gpuErrchk(cudaMalloc(&d_rate_matrix, 3 * t.getN()));

	// ----- allocate rates_constants arrays -----
	float * d_reaction_rates_constants;
	float * d_diffusion_rates_constants;
	gpuErrchk( cudaMalloc(&d_reaction_rates_constants, r.getR()) );
	gpuErrchk( cudaMalloc(&d_diffusion_rates_constants, r.getS()) );

	std::cout << " done!\n";

	// TODO: implementare update_rate_matrix su GPU

	gpuErrchk(cudaDeviceSynchronize());

}

}
