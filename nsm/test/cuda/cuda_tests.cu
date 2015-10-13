#include "cuda_tests.cuh"

void run_all()
{
	printf("---------- start CUDA tests ----------\n\n");

	printf("--- test_react_rates ...\n");
	float * d_react_rates_array;
	gpuErrchk(cudaMalloc(&d_react_rates_array, 3 * 4 * sizeof(float)));
	test_react_rates<<<4, 1>>>(d_react_rates_array);

	// TODO: we need to zero the memory between runs

	float want[3 * 4] =
		{ 16.0, 3.0, 2.0, 4.0, 60.0, 3.0, 3.0, 14.0, 128.0, 6.0, 8.0, 32.0 };
	float * got = new float[3 * 4];

	gpuErrchk(cudaMemcpy(got, d_react_rates_array, 3 * 4 * sizeof(float), cudaMemcpyDeviceToHost));

	for (int i = 0; i < 3 * 4; i++) {
		if (abs(got[i] - want[i]) > 1e-6) {
			printf("----- Failure in test_react_rates() -----\n");
			printf("position %d: got %.3f, want %.3f\n", i, got[i], want[i]);
			printf("-----------------------------------------\n\n");
		}
	}
	cudaDeviceReset();
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
		{ 0.5, 0.5, 0.5 };

	react_rates(state, reactants, sbc, spc, rc, rrc, react_rates_array);

}

__global__ void test_diff_rates(float * result)
{
	printf("--- test_diff_rates ...\n");

	// 0: 4 8 16 32
	// 1: 1 10 10 1
	// 2: 0 0 100 0
	int state[] =
		{ 4, 1, 0, 8, 10, 0, 16, 10, 100, 32, 1, 0 };
	int sbc = 3;
	int spc = 4;
	int sbi = 0;
	float drc[] =
		{ 0.1, 0.5, 0.5, 1.5 };

	diff_rates(state, sbc, spc, sbi, drc, result);
	float want[4] =
		{ 0.4, 4.0, 8.0, 48 };

	for (int i = 0; i < spc; i++) {
		if (abs(result[i] - want[i]) > 1e-6) {
			printf("----- Failure in test_diff_rates() -----\n");
			printf("specie %d: got %.3f, want %.3f\n", i, result[i], want[i]);
			printf("-----------------------------------------\n\n");
		}
	}

	sbi = 2;
	diff_rates(state, sbc, spc, sbi, drc, result);
	float want2[4] =
		{ 0, 0, 50.0, 0 };

	for (int i = 0; i < spc; i++) {
		if (abs(result[i] - want2[i]) > 1e-6) {
			printf("----- Failure in test_diff_rates() -----\n");
			printf("specie %d: got %.3f, want %.3f\n", i, result[i], want2[i]);
			printf("-----------------------------------------\n\n");
		}
	}
}

