#include <fstream>
#include <iostream>
#include <map>

#include "../include/cpp/Configuration.hpp"
#include "../include/cuda/driver.cuh"

int main(int argc, char * argv[])
{
	if (argc < 2) {
		std::cout << "Usage: ./pliscuda <configuration file>\n";
		return 1;
	} 

	std::cout << "\n    **   PLISCuda  **    \n\n";
	

	// check for CUDA enabled Device and print some info
	int nd;
	cudaGetDeviceCount(&nd);
	if(nd == 0 || nd > 99) { // sanity check 
		std::cout << "No CUDA-capable device detected.\n";
		return 1;
	}

	cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
	std::cout << "-- CUDA-capable device detected -- \n\n";
	std::cout << "  Device name: " << prop.name << "\n";
	std::cout << "  Total global memory: " << prop.totalGlobalMem / (1024*1024) << " MB\n";
	std::cout << "  Total shared memory per block: " << prop.sharedMemPerBlock / (1024)<< " KB\n";
	std::cout << "  Total registers per block: " <<  prop.regsPerBlock << "\n";
	std::cout << "  Number of multiprocessors: " << prop.multiProcessorCount << "\n\n";


	// open and parse configuration file
	std::cout << "-- Parsing input files -- \n\n";
	std::ifstream cf(argv[1]);
	if(!cf.good()) {
		std::cerr << "Couldn't load configuration file\n";
		return 1;
	}
	
	try {
		std::cerr << "  processing configuration file...\n";
		PLISCuda::Configuration conf(cf);

#ifdef DEBUG
		std::cout << conf << "\n";
#endif

		// start simulation
		PLISCuda::run_simulation(
			conf.getTopology(),
			conf.getState(),
			conf.getReactions(),
			conf.getReactionsConstantsArray(),
			conf.getDiffusionConstantsArray(),
			conf.getEndTime(),
			conf.getCompartments(),
			conf.getCompArray(),
			conf.getLogFreq()
			);
		
	} catch (std::invalid_argument &e) {
		std::cerr << "  parsing of configuration file failed with msg:\n";
		std::cerr << "\t" << e.what() << "\n";
	}

}

