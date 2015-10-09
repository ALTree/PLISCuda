#include "cuda_tests.cuh"

void run_all()
{
	driver<<<1, 1, 10>>>();
	cudaDeviceSynchronize();
}

__global__ void driver()
{
	test_react_rates();
}

__device__ void test_react_rates()
{
	// 0 0 1 ->
	// 0 2 0 ->
	// 1 0 1 ->
	int reactants[] = {
			0, 0, 1,
			0, 2, 0,
			1, 0, 1
	};
	int rc = 3;

	int state[] = {4, 8, 16};
	int sc = 1;
	int spc = 3;

	double rrc[] = { 0.2, 0.4, 0.5 };

	float * got = react_rates(reactants, rc, state, sc, spc, 0, rrc);
	float want[3] = {1.0, 11.200, 32.0};
	for(int i = 0; i < rc; i++) {
		if(abs(got[i] - want[i]) > 1e-6) {
			printf("----- Failure in test_react_rates() -----\n");
			printf("reaction %d: got %.3f, want %.3f\n", i, got[i], want[i]);
			printf("-----------------------------------------\n\n");
		}
	}
}
