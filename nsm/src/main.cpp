#include <fstream>
#include <iostream>

#include "../include/cpp/State.hpp"
#include "../include/cpp/Reactions.hpp"
#include "../include/cpp/Topology.hpp"

#include "../include/cpp/validation.hpp"
#include "../include/cuda/driver.cuh"

#include "../include/cuda/constants.cuh"

int main(int argc, char * argv[])
{
	if (argc < 8) {
		std::cout
				<< "Usage: ./nsm topology_file state_file reactions_file steps to_log subv_constants_file [constants_file]+ \n";
		return 1;
	}

#if LOG
	std::cout << "\n   ***   Start setup log   ***   \n\n";
	std::cout << "--- Parsing input files... ";
#endif

	// ---------- open and parse topology_file ----------
	std::ifstream topology_file(argv[1]);
	if (topology_file.bad()) {
		std::cerr << "couldn't load topology_file\n";
		return 1;
	}

	NSMCuda::Topology topology = NSMCuda::Topology();
	try {
		topology_file >> topology;
	} catch (std::invalid_argument &e) {
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
	} catch (std::invalid_argument &e) {
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
	} catch (std::invalid_argument &e) {
		std::cerr << "Parsing of reactions_file failed with\n";
		std::cerr << e.what() << "\n";
		return 1;
	}

	// --------- parse number of steps value ----------
	int steps = std::stoi(argv[4]);

	// --------- parse to_log file ----------
	std::ifstream log_file(argv[5]);
	bool * subv_to_log = new bool[topology.getN()];
	bool * spc_to_log = new bool[reactions.getS()];
	NSMCuda::read_log_data(log_file, subv_to_log, spc_to_log);

	// --------- parse subv <-> constants_set file ----------
	int * subv_constants = new int[topology.getN()];
	std::ifstream subv_file(argv[6]);
	NSMCuda::read_subv_constants(subv_file, subv_constants, topology.getN());

	// ---------- open and parse rate constants file(s) ----------
	int constants_files_count = argc - 7;

	float * reaction_rates_constants = new float[reactions.getR() * constants_files_count];
	float * diffusion_rates_constants = new float[reactions.getS() * constants_files_count];

	for (int i = 0; i < constants_files_count; i++) {
		std::ifstream constants_file(argv[i + 7]);
		NSMCuda::read_rates_constants(constants_file, &reaction_rates_constants[i * reactions.getR()],
				&diffusion_rates_constants[i * reactions.getS()], reactions.getR(), reactions.getS());
	}

#if LOG
	std::cout << " done!\n";
#endif

	NSMCuda::is_consistent(topology, initial_state, reactions);

#if LOG
	std::cout << "--- Consistency check... done!" << "\n\n";
#endif

	NSMCuda::run_simulation(topology, initial_state, reactions, reaction_rates_constants, diffusion_rates_constants,
			steps, constants_files_count, subv_constants);

}

