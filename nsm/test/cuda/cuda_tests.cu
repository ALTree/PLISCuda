#include "cuda_tests.cuh"

void run_all()
{
	printf("---------- start CUDA tests ----------\n\n");

	printf("--- test_react_rates ...\n");
	float * d_react_rates_array;
	gpuErrchk(cudaMalloc(&d_react_rates_array, 3 * 4 * sizeof(float)));
	gpuErrchk(cudaMemset(d_react_rates_array, 0, 3 * 4 * sizeof(float)));
	test_react_rates<<<4, 1>>>(d_react_rates_array);

	float want[3 * 4] =
		{ 16, 3, 2, 4, 180, 9, 9, 42, 512, 24, 32, 128 };
	float * got = new float[3 * 4];

	gpuErrchk(cudaMemcpy(got, d_react_rates_array, 3 * 4 * sizeof(float), cudaMemcpyDeviceToHost));

	for (int i = 0; i < 3 * 4; i++) {
		if (abs(got[i] - want[i]) > 1e-6) {
			printf("----- Failure in test_react_rates() -----\n");
			printf("position %d: got %.3f, want %.3f\n", i, got[i], want[i]);
			printf("-----------------------------------------\n\n");
		}
	}

	printf("--- test_diff_rates ...\n");
	float * d_diff_rates_array;
	gpuErrchk(cudaMalloc(&d_diff_rates_array, 3 * 4 * sizeof(float)));
	gpuErrchk(cudaMemset(d_diff_rates_array, 0, 3 * 4 * sizeof(float)));
	test_diff_rates<<<4, 1>>>(d_diff_rates_array);

	float want2[3 * 4] =
		{ 4, 1, 2, 4, 24, 6, 6, 12, 64, 12, 8, 16 };
	got = new float[3 * 4];

	gpuErrchk(cudaMemcpy(got, d_diff_rates_array, 3 * 4 * sizeof(float), cudaMemcpyDeviceToHost));

	for (int i = 0; i < 3 * 4; i++) {
		if (abs(got[i] - want2[i]) > 1e-6) {
			printf("----- Failure in test_diff_rates() -----\n");
			printf("position %d: got %.3f, want %.3f\n", i, got[i], want2[i]);
			printf("-----------------------------------------\n\n");
		}
	}


	cudaDeviceSynchronize();

	printf("\n----------  end  CUDA tests ----------\n");
}

__global__ void test_react_rates(float * react_rates_array)
{
	// 0: 8, 16, 32
	// 1: 2, 4, 6
	// 2: 4, 4, 4
	// 3: 8, 8, 8
	int state[] =
		{ 8, 2, 4, 8, 16, 4, 4, 8, 32, 6, 4, 8 };

	// 0 0 1 ->
	// 0 2 0 ->
	// 1 0 1 ->
	int reactants[] =
		{ 0, 0, 1, 0, 2, 0, 1, 0, 1 };

	int sbc = 4;
	int spc = 3;
	int rc = 3;

	float rrc[] =
		{ 0.5, 1.5, 2.0 };

	react_rates(state, reactants, sbc, spc, rc, rrc, react_rates_array);
}

__global__ void test_diff_rates(float * diff_rates_array)
{
	// 0: 8, 16, 32
	// 1: 2, 4, 6
	// 2: 4, 4, 4
	// 3: 8, 8, 8
	int state[] =
		{ 8, 2, 4, 8, 16, 4, 4, 8, 32, 6, 4, 8 };

	int sbc = 4;
	int spc = 3;

	float drc[] =
		{ 0.5, 1.5, 2.0 };

	diff_rates(state, sbc, spc, drc, diff_rates_array);
}

