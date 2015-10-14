#include <cuda_runtime.h>

#include "../../include/cuda/nsm_driver.hpp"
#include "../../include/cuda/rates.cuh"
#include "../../include/cuda/cuda_utils.cuh"

namespace NSMCuda {

void nsm(Topology t, State s, Reactions r, float * h_rrc, float * h_drc)
{

	int sbc = t.getN();
	int spc = s.getS();
	int rc = r.getR();

	std::cout << "----- Allocating GPU memory ...";

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
	gpuErrchk(cudaMalloc(&d_reactants, rc * sizeof(int)));
	gpuErrchk(cudaMalloc(&d_products, rc * sizeof(int)));
	gpuErrchk(cudaMemcpy(d_reactants, h_reactants, rc * sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_products, h_products, rc * sizeof(int), cudaMemcpyHostToDevice));

	// ----- allocate and memcpy topology array -----
	int * h_topology = t.getNeighboursArray();

	int * d_topology;
	gpuErrchk(cudaMalloc(&d_topology, 6 * sbc * sizeof(int)));
	gpuErrchk(cudaMemcpy(d_topology, h_topology, 6 * sbc * sizeof(int), cudaMemcpyHostToDevice));

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
	gpuErrchk(cudaMalloc(&d_react_rates_array, sbc * rc));
	gpuErrchk(cudaMalloc(&d_diff_rates_array, sbc * spc));

	// zero GPU memory, just to be sure
	// TODO: remove(?)
	gpuErrchk(cudaMemset(d_rate_matrix, 0, 3 * sbc));
	gpuErrchk(cudaMemset(d_react_rates_array, 0, sbc * rc));
	gpuErrchk(cudaMemset(d_diff_rates_array, 0, sbc * spc));

	std::cout << " done!\n";

	gpuErrchk(cudaDeviceSynchronize());

}

}
