#ifndef CPP_UTILS_HPP_
#define CPP_UTILS_HPP_

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

// this NEEDS to stay here: it's to be used in host code
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort)
			exit(code);
	}
}

#endif /* CPP_UTILS_HPP_ */
