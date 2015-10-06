#include <cassert>

#include "../include/State.hpp"
#include "../include/Reactions.hpp"
#include "../include/Topology.hpp"

namespace NSMCuda {

bool are_consistent(Topology t, State s, Reactions r)
{
	// checks on subvolumes number
	assert(t.getN() == s.getN());
	return true;

}


}
