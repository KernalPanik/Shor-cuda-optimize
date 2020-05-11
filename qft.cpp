#include "qft.h"

#include <stdio.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <device_launch_parameters.h>
#include <ostream>
#include <iostream>
#include <cuComplex.h>
#include <complex>

#include <chrono>
#include "types.h"
#include "quantumRegister.h"

#define RANK 2

void CuFFT_QFT_Calc(amp *i, amp *o)
{

}


state_map CalculateQFT(std::vector<string> statenames, std::vector<amp> statevals)
{
	cufftHandle plan;
	cufftComplex* inputData;
	cufftComplex* outputData;
	std::vector<std::complex<double>> fcompl;


	for (auto v : statevals)
	{
		fcompl.push_back(v);
	}

	if (cufftPlan1d(&plan, statevals.size(), CUFFT_C2C, 1) != CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT error: Plan creation failed");
	}


	cudaMallocManaged((void**)&inputData, sizeof(cufftComplex) * statevals.size() * 1);
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Cuda error: Failed to allocate memory for input data\n");
	}

	int i = 0;
	for (auto v : fcompl)
	{
		inputData[i].x = v.real();
		inputData[i].y = v.imag();
 		i++;
	}

	cudaMallocManaged((void**)&outputData, sizeof(cufftComplex) * statevals.size() * 1);
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Cuda error: Failed to allocate memory for output data\n");
	}

	if (cufftExecC2C(plan, inputData, outputData, CUFFT_FORWARD) != CUFFT_SUCCESS)
	{
		fprintf(stderr, "CUFFT error: ExecC2C failed");
	}

	if (cudaDeviceSynchronize() != cudaSuccess) {
		fprintf(stderr, "Cuda error: Failed to synchronize\n");
	}

	cufftDestroy(plan);
	cudaFree(inputData);

	std::map<std::string, amp> a;

	for (int i = 0; i < statevals.size(); i++)
	{
		amp temp(outputData[i].x, outputData[i].y);

		a.insert({ statenames[i], temp });
	}
	cudaFree(outputData);

	return a;
}
