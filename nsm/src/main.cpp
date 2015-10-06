#include <fstream>
#include <iostream>

#include "../include/State.hpp"
#include "../include/Reactions.hpp"
#include "../include/Topology.hpp"


int main(int argc, char ** argv)
{
	if (argc < 4) {
		std::cout << "Usage: ./nsm topology_file state_file reactions_file\n";
		return 1;
	}

	// ---------- open and parse topology_file ----------

	std::ifstream topology_file(argv[1]);
	if (topology_file.bad()) {
		std::cerr << "couldn't load topology_file\n";
		return 1;
	}

	NSMCuda::Topology topology = NSMCuda::Topology();
	try {
		topology_file >> topology;
	} catch(std::invalid_argument &e) {
		std::cerr << "Parsing of topology_file failed with\n";
		std::cerr << e.what() << "\n";
		return 1;
	}


	// ---------- open and parse state_file ----------

	std::ifstream state_file(argv[2]);
	if (state_file.bad()) {
		std::cerr << "couldn't load state_file\n";
		return 1;
	}

	NSMCuda::State initial_state = NSMCuda::State();
	try {
		state_file >> initial_state;
	} catch(std::invalid_argument &e) {
		std::cerr << "Parsing of state_file failed with\n";
		std::cerr << e.what() << "\n";
		return 1;
	}


	// ---------- open and parse reactions_file ----------

	std::ifstream reactions_file(argv[3]);
	if (reactions_file.bad()) {
		std::cerr << "couldn't load reactions_file\n";
		return 1;
	}

	NSMCuda::Reactions reactions = NSMCuda::Reactions();
	try {
		reactions_file >> reactions;
	} catch(std::invalid_argument &e) {
		std::cerr << "Parsing of reactions_file failed with\n";
		std::cerr << e.what() << "\n";
		return 1;
	}


	// TODO: implement sanity checks on topology, initial_state and reactions

}

