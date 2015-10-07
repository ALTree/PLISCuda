#include <cuda_runtime.h>

#include "../include/nsm_driver.hpp"
#include "../include/nsm.cuh"

namespace NSMCuda {

void nsm(Topology t, State s, Reactions r)
{

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
	int * h_offsets_array = t.getOffsetArray();
	int topology_len = t.getNeighboursLength();
	int offsets_len = t.getN();

	int * d_topology_array, *d_offsets_array;
	gpuErrchk(cudaMalloc(&d_topology_array, topology_len));
	gpuErrchk(cudaMalloc(&d_offsets_array, offsets_len));

	// ----- allocate rate matrix -----
	float * d_rate_matrix;
	gpuErrchk(cudaMalloc(&d_rate_matrix, 3 * t.getN()));

	std::cout << " done!\n";

	foo();

	gpuErrchk( cudaDeviceSynchronize() );



}

}
