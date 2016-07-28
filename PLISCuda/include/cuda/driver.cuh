#ifndef NSM_DRIVER_HPP_
#define NSM_DRIVER_HPP_

#include <ctime>
#include <fstream>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/logical.h>
#include <curand_kernel.h>
#include <iomanip>

#include "../cpp/State.hpp"
#include "../cpp/Reactions.hpp"
#include "../cpp/Topology.hpp"

#include "../cpp/cpp_utils.hpp"

#include "cuda_utils.cuh"
#include "rates.cuh"
#include "ssa.cuh"
#include "leap.cuh"
#include "constants.cuh"
#include "log.cuh"

struct ToLog {
  unsigned int * subv;
  int subv_len;
  bool * spc;
  int spc_len;
  float freq;
};

namespace PLISCuda {

  // Run a simulation with the given topology, initial state, set of reactions, and
  // react and diffusion rates constants, until simulation time exceeds *endTime*.
  void run_simulation(Topology t, State s, Reactions r, float * h_rrc, float * h_drc, float endTime,
					  int constants_files_count, int * subv_constants, float logFreq);

  // get index of min value in tau array
  int h_get_min_tau(thrust::device_vector<float> &tau);

  // Print utils for logging
  void print_state(int * h_state, int spc, int sbc, float current_time);
  void print_rate_matrix(float * h_rate_matrix, int sbc);
  void print_tau(thrust::device_vector<float>, int sbc);
  void print_leap_array(char * d_leap, int sbc);
  std::string print_state_snapshot(int * h_state, int spc, int sbc, float current_time);

}

#endif /* NSM_DRIVER_HPP_ */
