#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <iostream>

using namespace std;

class Base
{
public:
    __host__ __device__ Base() {}
    __host__ __device__ ~Base() {}

    __host__ __device__ virtual bool fun(int index)
    {
        printf("No.%d Base.fun(int)\n", index);
        return true;
    }

    int value;
};

class Derived :public Base
{
public:
    __host__ __device__ Derived() {}
    __host__ __device__ ~Derived() {}
    __host__ __device__ virtual bool fun(int index)
    {
        printf("No.%d Derived.fun(int)\n", index);
        return true;
    }
};

__global__ void TestFun(bool *d_result, Base **d_derives, size_t size)
{
    int myId = blockDim.x*blockIdx.x + threadIdx.x;
    if (myId >= size)
        return;

    (*d_derives)[myId].value = myId;
    d_result[myId] = (*d_derives)[myId].fun(myId);
}

// 在device上申请Derived类的实例
__global__ void AllocateOnDevice(Base **d_b, const size_t size)
{
    *d_b = new Derived[size];
}

// 释放之前在device上申请Derived的实例
__global__ void DeleteOnDevice(Base **d_b)
{
    delete[] (*d_b);
    d_b = nullptr;
}

int main()
{
    bool *d_result = nullptr;
    cudaMalloc((void **)&d_result, sizeof(bool) * 10);

    Base **d_derives = nullptr;
    cudaMalloc((void **)&d_derives, sizeof(Base *));
    AllocateOnDevice <<<1, 1 >>> (d_derives, 254000);       // 在device上申请Derived类的实例

    TestFun << <1, 10 >> > (d_result, d_derives, 10);

    cudaDeviceSynchronize();
    cudaGetLastError();

    cudaFree(d_result);
    DeleteOnDevice <<<1, 1 >>>(d_derives);              // 释放之前在device上申请Derived的实例

    return 0;
}