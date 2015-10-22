#include "cuda_tests.cuh"

int main()
{
	printf("---------- start CUDA tests ----------\n\n");
	run_rates_tests();
	run_rand_tests();
	printf("\n----------  end  CUDA tests ----------\n");
}
