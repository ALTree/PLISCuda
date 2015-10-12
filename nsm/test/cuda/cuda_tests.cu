#include "cuda_tests.cuh"

void run_all()
{
	printf("---------- start CUDA tests ----------\n\n");

	test_react_rates<<<1, 1, 3 * sizeof(float)>>>();
	test_diff_rates<<<1, 1, 4 * sizeof(float)>>>();
	cudaDeviceSynchronize();

	printf("\n----------  end  CUDA tests ----------\n");
}

__global__ void test_react_rates()
{
	printf("--- test_react_rates ...\n");
	// 0 0 1 ->
	// 0 2 0 ->
	// 1 0 1 ->
	int reactants[] =
		{ 0, 0, 1, 0, 2, 0, 1, 0, 1 };
	int rc = 3;

	int state[] =
		{ 4, 8, 16 };
	int sc = 1;
	int spc = 3;

	float rrc[] =
		{ 0.2, 0.4, 0.5 };

	float * got = react_rates(reactants, rc, state, sc, spc, 0, rrc);
	float want[3] =
		{ 3.2, 11.2, 32.0 };
	for (int i = 0; i < rc; i++) {
		if (abs(got[i] - want[i]) > 1e-6) {
			printf("----- Failure in test_react_rates() -----\n");
			printf("reaction %d: got %.3f, want %.3f\n", i, got[i], want[i]);
			printf("-----------------------------------------\n\n");
		}
	}
}

__global__ void test_diff_rates()
{
	printf("--- test_diff_rates ...\n");

	// 0: 4 8 16 32
	// 1: 1 10 10 1
	// 2: 0 0 100 0
	int state[] =
		{ 4, 1, 0, 8, 10, 0, 16, 10, 100, 32, 1, 0 };
	int sc = 3;
	int spc = 4;

	float drc[] =
		{ 0.1, 0.5, 0.5, 1.5 };

	float * got = diff_rates(state, sc, spc, 0, drc);
	float want[4] =
		{ 0.4, 4.0, 8.0, 48 };

	for (int i = 0; i < spc; i++) {
		if (abs(got[i] - want[i]) > 1e-6) {
			printf("----- Failure in test_diff_rates() -----\n");
			printf("specie %d: got %.3f, want %.3f\n", i, got[i], want[i]);
			printf("-----------------------------------------\n\n");
		}
	}

	got = diff_rates(state, sc, spc, 2, drc);
	float want2[4] =
		{ 0, 0, 50.0, 0 };

	for (int i = 0; i < spc; i++) {
		if (abs(got[i] - want2[i]) > 1e-6) {
			printf("----- Failure in test_diff_rates() -----\n");
			printf("specie %d: got %.3f, want %.3f\n", i, got[i], want2[i]);
			printf("-----------------------------------------\n\n");
		}
	}

}

