#ifndef NSM_DRIVER_HPP_
#define NSM_DRIVER_HPP_

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/logical.h>
#include <curand_kernel.h>

#include "../cpp/State.hpp"
#include "../cpp/Reactions.hpp"
#include "../cpp/Topology.hpp"

#include "../cpp/cpp_utils.hpp"

#include "cuda_utils.cuh"
#include "rates.cuh"
#include "nsm.cuh"
#include "leap.cuh"
#include "constants.cuh"

namespace NSMCuda {

// Run a simulation with the given topology, initial state, set of reactions, and
// react and diffusion rates constants, for *steps* steps.
void run_simulation(Topology t, State s, Reactions r, float * h_rrc, float * h_drc, int steps);

// Print utils for logging
void print_state(int * h_state, int spc, int sbc);
void print_rate_matrix(float * h_rate_matrix, int sbc);
void print_tau(thrust::device_vector<float>, int sbc);

}

#endif /* NSM_DRIVER_HPP_ */
