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

}

#endif /* NSM_DRIVER_HPP_ */
