#include "cuda_tests.cuh"

void run_all()
{
	printf("---------- start CUDA tests ----------\n\n");

	// ----- test_react_rates -----
	float * react_rates_result;
	cudaMalloc(&react_rates_result, 3 * sizeof(float));
	test_react_rates<<<1, 1>>>(react_rates_result);

	// ----- test_diff_rates -----
	float * diff_rates_result;
	cudaMalloc(&diff_rates_result, 4 * sizeof(float));
	test_diff_rates<<<1, 1>>>(diff_rates_result);

	float * rate_matrix;
	cudaMalloc(&rate_matrix, 3 * 3);
	test_rate_matrix_row<<<1, 1, 32 * sizeof(float)>>>(rate_matrix);

	cudaDeviceSynchronize();

	printf("\n----------  end  CUDA tests ----------\n");
}

__global__ void test_react_rates(float * result)
{
	printf("--- test_react_rates ...\n");

	int state[] =
		{ 4, 8, 16 };

	// 0 0 1 ->
	// 0 2 0 ->
	// 1 0 1 ->
	int reactants[] =
		{ 0, 0, 1, 0, 2, 0, 1, 0, 1 };

	int sbc = 1;
	int spc = 3;
	int rc = 3;
	int sbi = 0;

	float rrc[] =
		{ 0.2, 0.4, 0.5 };

	react_rates(state, reactants, sbc, spc, rc, sbi, rrc, result);
	float want[3] =
		{ 3.2, 11.2, 32.0 };
	for (int i = 0; i < rc; i++) {
		if (abs(result[i] - want[i]) > 1e-6) {
			printf("----- Failure in test_react_rates() -----\n");
			printf("reaction %d: got %.3f, want %.3f\n", i, result[i], want[i]);
			printf("-----------------------------------------\n\n");
		}
	}
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

__global__ void test_rate_matrix_row(float * rate_matrix)
{
	/*

	 // 0: 2 2 2
	 // 1: 4 4 8
	 // 2: 8 8 4
	 int state[9] =
	 { 2, 4, 8, 2, 4, 8, 2, 8, 4, };

	 // 0 0 1 ->
	 // 0 2 0 ->
	 // 1 0 1 ->
	 int reactants[] =
	 { 0, 0, 1, 0, 2, 0, 1, 0, 1 };

	 int sbc = 3;
	 int spc = 3;
	 int rc = 3;

	 float rrc[3] =
	 { 1, 2, 3 };

	 float drc[3] =
	 { 0.5, 0.5, 0.5 };

	 rate_matrix_row(state, reactants, sbc, spc, rc, 0, rate_matrix, rrc, drc);

	 float want[3] =
	 { 1, 2, 3 };

	 for (int i = 0; i < 3; i++) {
	 if (abs(rate_matrix[sbc * i] - want[i]) > 1e-6) {
	 printf("----- Failure in test_rate_matrix_row() -----\n");
	 printf("column %d: got %.3f, want %.3f\n", i, rate_matrix[sbc * i], want[i]);
	 printf("-----------------------------------------\n\n");
	 }
	 }
	 */
}

