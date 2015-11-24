#ifndef NSM_DRIVER_HPP_
#define NSM_DRIVER_HPP_

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
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

void nsm(Topology, State, Reactions, float *, float *);

void print_state(int * h_state, int spc, int sbc);
void print_rate_matrix(float * h_rate_matrix, int sbc);
void print_tau(thrust::device_vector<float>, int sbc);

}

#endif /* NSM_DRIVER_HPP_ */
