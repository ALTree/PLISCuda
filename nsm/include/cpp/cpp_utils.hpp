#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <string>
#include <sstream>
#include <vector>

#include <cuda_runtime.h>

namespace NSMCuda {

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);

std::vector<std::string> split(const std::string &s, char delim);

bool read_rates_constants(std::istream& is, float * reaction_rates_constants, float * diffusion_rates_constants,
		int reactions_count, int species_count);

}

#endif /* UTILS_HPP_ */
