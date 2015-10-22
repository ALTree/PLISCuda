#include "cuda_tests.cuh"

void run_rates_tests()
{

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

	printf("--- test_update_rate_matrix ...\n");
	float * d_rate_matrix;

	gpuErrchk(cudaMalloc(&d_rate_matrix, 3 * 3 * sizeof(float)));
	gpuErrchk(cudaMalloc(&d_react_rates_array, 3 * 2 * sizeof(float)));
	gpuErrchk(cudaMalloc(&d_diff_rates_array, 3 * 3 * sizeof(float)));

	gpuErrchk(cudaMemset(d_rate_matrix, 0, 3 * 3 * sizeof(float)));
	gpuErrchk(cudaMemset(d_react_rates_array, 0, 3 * 2 * sizeof(float)));
	gpuErrchk(cudaMemset(d_diff_rates_array, 0, 3 * 3 * sizeof(float)));

	test_update_rate_matrix<<<3, 1>>>(d_rate_matrix, d_react_rates_array, d_diff_rates_array);

	float want3[] =
		{ 60, 14, 60, 48, 48, 48, 108, 62, 108 };

	gpuErrchk(cudaMemcpy(got, d_rate_matrix, 3 * 3 * sizeof(float), cudaMemcpyDeviceToHost));

	for (int i = 0; i < 3 * 3; i++) {
		if (abs(got[i] - want3[i]) > 1e-6) {
			printf("----- Failure in test_update_rate_matrix() -----\n");
			printf("position %d: got %.3f, want %.3f\n", i, got[i], want3[i]);
			printf("-----------------------------------------\n\n");
		}
	}

	cudaDeviceSynchronize();
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

__global__ void test_update_rate_matrix(float * rate_matrix, float * react_rates_array, float * diff_rates_array)
{
	// 0: 1
	// 1: 0 2
	// 2: 1
	int topology[] =
		{ 1, -1, -1, -1, -1, -1, 0, 2, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1 };

	// 0: 8 8 8
	// 1: 4 4 4
	// 2: 8 8 8
	int state[] =
		{ 8, 4, 8, 8, 4, 8, 8, 4, 8 };

	// 1 0 1 ->
	// 0 2 0 ->
	int reactants[] =
		{ 1, 0, 0, 2, 1, 0 };

	int sbc = 3;
	int spc = 3;
	int rc = 2;

	float rrc[] =
		{ 0.5, 1 };

	float drc[] =
		{ 2, 2, 2 };

	react_rates(state, reactants, sbc, spc, rc, rrc, react_rates_array);
	diff_rates(state, sbc, spc, drc, diff_rates_array);
	update_rate_matrix(topology, sbc, spc, rc, rate_matrix, react_rates_array, diff_rates_array);
}

// --------------------------------------------------

void run_rand_tests()
{
	printf("--- test_fill_tau_array ...\n");

	int sbc = 32;
	thrust::device_vector<float> tau(sbc);
	fill_tau_array<<<1, sbc>>>(thrust::raw_pointer_cast(tau.data()), tau.size());

	for (int i = 0; i < tau.size(); i++) {
		if (tau[i] == 0) {
			printf("----- Failure in test_fill_tau_array() -----\n");
			printf("position %d: found a zero tau\n", i);
			printf("-----------------------------------------\n\n");
		}
	}

	printf("--- test_choose_random_reaction ...\n");
	test_choose_random_reaction<<<1, 4>>>();


}

__global__ void test_choose_random_reaction()
{
	int sbc = 4;
	int rc = 3;

	float rate_matrix[] = {
			1.7, 2.3, 3.2, 4.3, // sums of reaction rates (R)
			1.3, 2.7, 3.8, 4.7, // sums of diffusion rates (D)
			3.0, 5.0, 7.0, 9.0  // R + D
	};

	float react_rates_array[] = {
			1.0, 1.3, 1.1, 0.3,       // reaction 1
			0.4, 0.4, 1.1, 2.3,       // reaction 2
			0.3, 0.6, 1.0, 1.7        // reaction 3
	};

	int sbi = blockIdx.x * blockDim.x + threadIdx.x;

	curandState s;
	curand_init(sbi, 0, 0, &s);
	for(int i = 0; i < 32; i++) {
		int ri = choose_rand_reaction(sbc, rc, rate_matrix, react_rates_array, curand_uniform(&s));
		if(ri < 0 || ri >= rc) {
			printf("----- Failure in test_fill_tau_array() -----\n");
			printf("sbi = %d, random reaction index is off: %d\n", i, ri);
			printf("-----------------------------------------\n\n");
		}
	}

}
