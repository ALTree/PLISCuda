#include "../include/validation.hpp"

namespace NSMCuda {

void is_consistent(NSMCuda::Topology t, NSMCuda::State s, NSMCuda::Reactions r)
{
	// checks on subvolumes number
	assert(t.getN() == s.getN());

	// check on species number
	assert(s.getS() == r.getS());

}

}
