#include <iostream>
#include <cuda_runtime.h>

// Declare a __device__ variable
__device__ int deviceVar[3];

__global__ void kernel()
{
    // Modify the __device__ variable
    deviceVar[0] = 1;
    deviceVar[1] = 2;
    deviceVar[2] = 3;
}

int main()
{
    // Launch the kernel
    kernel<<<1, 1>>>();

    // Copy the __device__ variable to host
    int hostVar[3];
    cudaMemcpyFromSymbol(hostVar, deviceVar, 3 * sizeof(int), 0, cudaMemcpyDeviceToHost);

    // Print the value of the __device__ variable
    std::cout << "Values of deviceVar: " << hostVar[0] << ", " << hostVar[1] << ", " << hostVar[2] << std::endl;

    return 0;
}