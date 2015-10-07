#ifndef NSM_DRIVER_HPP_
#define NSM_DRIVER_HPP_

#include "State.hpp"
#include "Reactions.hpp"
#include "Topology.hpp"

#include "utils.hpp"

namespace NSMCuda {

void nsm(Topology t, State s, Reactions r);

}

#endif /* NSM_DRIVER_HPP_ */
