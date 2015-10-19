#ifndef NSM_DRIVER_HPP_
#define NSM_DRIVER_HPP_

#include "../cpp/State.hpp"
#include "../cpp/Reactions.hpp"
#include "../cpp/Topology.hpp"

#include "../cpp/cpp_utils.hpp"

namespace NSMCuda {

void nsm(Topology, State, Reactions, float *, float *);

}

#endif /* NSM_DRIVER_HPP_ */
