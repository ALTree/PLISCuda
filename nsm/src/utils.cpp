#include "../include/utils.hpp"

namespace NSMCuda {

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

std::vector<std::string> split(const std::string &s, char delim)
{
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

bool read_rates_constants(std::istream& is, float * reaction_rates_constants, float * diffusion_rates_constants,
		int reactions_count, int species_count)
{
	std::string line;

	std::getline(is, line);
	std::vector<std::string> rate_constants = NSMCuda::split(line, ' ');
	for(int i = 0; i < reactions_count; i++) {
		reaction_rates_constants[i] = std::stof(rate_constants.at(i));
	}

	std::getline(is, line);
	rate_constants = NSMCuda::split(line, ' ');
	for(int i = 0; i < species_count; i++) {
		diffusion_rates_constants[i] = std::stof(rate_constants.at(i));
	}

	return true;


}

}
