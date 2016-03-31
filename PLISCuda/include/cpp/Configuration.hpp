#ifndef CONFIGURATION_HPP_
#define CONFIGURATION_HPP_

#include <map>

#include "Topology.hpp"

namespace NSMCuda {

class Configuration {

	Topology topology;
	Reactions reactions;
	State initial_state;
	int steps;
	std::vector<std::vector<float>> reactions_constants;
	std::vector<std::vector<float>> diffusion_constants;

public:

	// ---------- Constructor ----------
	Configuration(std::ifstream& config_file)
	{
		std::string line;
		std::ifstream fs;
		
		// parse topology line
		std::getline(config_file, line, '=');
		if (line != "topology") {
			throw std::invalid_argument("config file: need a topology file");
		}
		std::getline(config_file, line);
		fs.open(line);
		if (!fs.good()) {
			throw std::invalid_argument("config file: " + line + " not found");
		}
		fs >> topology;
		fs.close(); fs.clear();

		// parse reactions line
		std::getline(config_file, line, '=');
		if (line != "reactions") {
			throw std::invalid_argument("config file: need a reactions file");
		}
		std::getline(config_file, line);
		fs.open(line);
		if (!fs.good()) {
			throw std::invalid_argument("config file: " + line + " not found");
		}
		fs >> reactions;
		
		
		
	}

	// ---------- Destructor ----------
	~Configuration()
	{
		
	}


	// ---------- Getters ----------
	Topology getTopology()
	{
		return topology;
	}

	Reactions getReactions()
	{
		return reactions;
	}

	// int getS() const
	// {
	// 	return s;
	// }

	// int * getReactants() const
	// {
	// 	return reactants;
	// }

	// int * getProducts() const
	// {
	// 	return products;
	// }

};

}

#endif /* CONFIGURATION_HPP_ */
