#ifndef CONFIGURATION_HPP_
#define CONFIGURATION_HPP_

#include <map>

#include "Topology.hpp"
#include "State.hpp"
#include "Reactions.hpp"

namespace PLISCuda {

class Configuration {

	Topology topology;
	Reactions reactions;
	State initial_state;
	float end_time;
	float log_freq;
	int compartments;
	int * comp_array;
	float * reactions_constants;
	float * diffusion_constants;

public:

	// ---------- Constructor ----------
	Configuration(std::ifstream& config_file)
	{
		std::string line;
		std::ifstream fs;
		
		// parse topology line
		std::getline(config_file, line, '=');
		if (line != "topology") {
			throw std::invalid_argument("need a topology file");
		}
		std::getline(config_file, line);
		fs.open(line);
		if (!fs.good()) {
			throw std::invalid_argument("file " + line + " not found");
		}
		fs >> topology;
		fs.close(); fs.clear();

		// parse reactions line
		std::getline(config_file, line, '=');
		if (line != "reactions") {
			throw std::invalid_argument("need a reactions file");
		}
		std::getline(config_file, line);
		fs.open(line);
		if (!fs.good()) {
			throw std::invalid_argument("file " + line + " not found");
		}
		fs >> reactions;
		fs.close(); fs.clear();

		// parse initial state line
		std::getline(config_file, line, '=');
		if (line != "initialState") {
			throw std::invalid_argument("need a state file");
		}
		std::getline(config_file, line);
		fs.open(line);
		if (!fs.good()) {
			throw std::invalid_argument("file " + line + " not found");
		}
		fs >> initial_state;
		fs.close(); fs.clear();
		
		// parse end time line
		std::getline(config_file, line, '=');
		if (line != "endTime") {
			throw std::invalid_argument("need an endtime value");
		}
		std::getline(config_file, line);
		end_time = std::stof(line);

		// parse log freq line
		std::getline(config_file, line, '=');
		if (line != "logFreq") {
			throw std::invalid_argument("need a logFreq value (0 = no logging)");
		}
		std::getline(config_file, line);
		log_freq = std::stof(line);

		// parse compartments line
		std::getline(config_file, line, '=');
		if (line != "compartments") {
			throw std::invalid_argument("need a compartments value");
		}
		std::getline(config_file, line);
		compartments = std::stoi(line);
		if(compartments < 1) {
			throw std::invalid_argument("compartments count < 1");
		}

		// if compartments > 1 we need a subvolume <-> compartments file
		if(compartments > 1) {
			std::getline(config_file, line, '=');
			if (line != "compartmentsFile") {
				throw std::invalid_argument("need a compartments file");
			}
			std::getline(config_file, line);
			fs.open(line);
			if (!fs.good()) {
				throw std::invalid_argument("file " + line + " not found");
			}

			int subvs = topology.getN();
			comp_array = new int[subvs];
			PLISCuda::read_subv_constants(fs, comp_array, subvs);
			fs.close(); fs.clear();
		} else {
			// generate faux compartmentsFile on the fly
			int subvs = topology.getN();
			comp_array = new int[subvs](); // zero-initialized
		}

		// parse constants file(s)
		int rcount = reactions.getR();
		int scount = reactions.getS();
		reactions_constants = new float[rcount * compartments];
		diffusion_constants = new float[scount * compartments];

		// we need a file for each compartment
		for(int i = 0; i < compartments; i++) {
			std::getline(config_file, line, '=');
			if (line != "constants") {
				throw std::invalid_argument("need a constants file (" + std::to_string(i+1) + ")");
			}
			std::getline(config_file, line);
			fs.open(line);
			if (!fs.good()) {
				throw std::invalid_argument("file " + line + " not found");
			}

			PLISCuda::read_rates_constants(fs,
										  &reactions_constants[i * rcount],
										  &diffusion_constants[i * scount],
										  rcount, scount);
			fs.close(); fs.clear();
		}
	}

	// ---------- Destructor ----------
	~Configuration()
	{
		delete[] comp_array;
		delete[] reactions_constants;
		delete[] diffusion_constants;
	}


	// ---------- Getters ----------
	Topology& getTopology()
	{
		return topology;
	}

	Reactions& getReactions()
	{
		return reactions;
	}

	State& getState()
	{
		return initial_state;
	}

	float getEndTime() const
	{
		return end_time;
	}

	float getLogFreq() const
	{
		return log_freq;
	}

	int getCompartments() const
	{
		return compartments;
	}

	int * getCompArray() const
	{
		return comp_array;
	}

	float * getReactionsConstantsArray() const
	{
		return reactions_constants;
	}

	float * getDiffusionConstantsArray() const
	{
		return diffusion_constants;
	}
	

};

inline std::ostream& operator<<(std::ostream& os, Configuration& c)
{
	os << "--- Configuration --- " << "\n\n";
	os << "  " << c.getTopology() << "\n";
	os << "  " << c.getReactions() << "\n";
	os << "  " << c.getState() << "\n";
	os << "  endTime = " << c.getEndTime() << "\n\n";
	os << "  logFreq = " << c.getLogFreq() << "\n\n";
	int comp = c.getCompartments();
	os << "  compartments = " << comp << "\n";
	if(c.getCompartments() > 1) {
		for (int i = 0; i < c.getTopology().getN(); i++) {
			os << c.getCompArray()[i] << " ";
		}
		os << "\n";
	}

	os << "\n  reactions constants:\n";
	os << "  ";
	for (int i = 0; i < comp * c.getReactions().getR(); i++) {
		os << c.getReactionsConstantsArray()[i] << " ";
	}

	os << "\n\n  diffusion constants:\n";
	os << "  ";
	for (int i = 0; i < comp * c.getReactions().getS(); i++) {
		os << c.getDiffusionConstantsArray()[i] << " ";
	}

	os << "\n";
	
	return os;
}

}

#endif /* CONFIGURATION_HPP_ */
