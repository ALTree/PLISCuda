#ifndef CPP_UTILS_HPP_
#define CPP_UTILS_HPP_

#include <string>
#include <sstream>
#include <vector>

#include <cuda_runtime.h>

namespace NSMCuda {

inline std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

inline std::vector<std::string> split(const std::string &s, char delim)
{
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

inline void read_rates_constants(std::istream& is, float * reaction_rates_constants, float * diffusion_rates_constants,
		int reactions_count, int species_count)
{
	std::string line;

	std::getline(is, line);
	std::vector<std::string> rate_constants = NSMCuda::split(line, ' ');
	for (int i = 0; i < reactions_count; i++) {
		reaction_rates_constants[i] = std::stof(rate_constants.at(i));
	}

	std::getline(is, line);
	rate_constants = NSMCuda::split(line, ' ');
	for (int i = 0; i < species_count; i++) {
		diffusion_rates_constants[i] = std::stof(rate_constants.at(i));
	}
}

inline void read_subv_constants(std::istream& is, int * subv_constants, uint sbc)
{
	std::string line;

	for (int i = 0; i < sbc; i++) {
		std::getline(is, line, ':');
		std::getline(is, line);
		subv_constants[i] = std::stoi(line);
	}
}

inline float read_log_data(std::istream& is, bool * log_subv, bool * log_specie)
{
	std::string line;

	// parse log frequency
	std::getline(is, line, ':');
	std::getline(is, line);
	float freq = std::stof(line);

	// parse subvolumes line
	std::getline(is, line, ' ');
	std::getline(is, line);
	std::vector<std::string> subv = NSMCuda::split(line, ' ');
	for (auto &i : subv) {
		log_subv[std::stoi(i)] = true;
	}

	// parse species line
	std::getline(is, line, ' ');
	std::getline(is, line);
	std::vector<std::string> spc = NSMCuda::split(line, ' ');
	for (auto &i : spc) {
		log_specie[std::stoi(i)] = true;
	}

	return freq;
}

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
