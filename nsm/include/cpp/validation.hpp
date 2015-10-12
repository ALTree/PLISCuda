#ifndef VALIDATION_HPP_
#define VALIDATION_HPP_

#include <cassert>

#include "State.hpp"
#include "Topology.hpp"
#include "Reactions.hpp"

namespace NSMCuda {

void is_consistent(NSMCuda::Topology t, NSMCuda::State s, NSMCuda::Reactions r);

}

#endif /* VALIDATION_HPP_ */
