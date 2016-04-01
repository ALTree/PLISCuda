#ifndef VALIDATION_HPP_
#define VALIDATION_HPP_

#include <cassert>

#include "State.hpp"
#include "Topology.hpp"
#include "Reactions.hpp"

namespace PLISCuda {

inline void is_consistent(PLISCuda::Topology t, PLISCuda::State s, PLISCuda::Reactions r)
{
	// checks on subvolumes number
	assert(t.getN() == s.getN());

	// check on species number
	assert(s.getS() == r.getS());
}

}

#endif /* VALIDATION_HPP_ */
