#include "../include/nsm.cuh"

__global__ void test()
{
	printf("hello!\n");
}

void foo()
{
	test<<<10,1>>>();
}
