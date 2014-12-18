#ifndef PETSC_WIND_CUH
#define PETSC_WIND_CUH

#include <cuda.h>
#include <cuda_runtime.h>

void initGPUWindResources(cudaArray **windVelArrayDevice, 
	int dimx, int dimy, int dimz, int tdim);
void windToGPU(float4 *windVel, cudaArray *windVelArrayDevice, 
	int dimx, int dimy, int dimz);
void freeGPUWindResources(cudaArray *array);
void getTerrainMapFromGPU(float4 *terrainMap, int tdim);

#endif