#include "PetscWind.cuh"
#include "CudaHelpers.cuh"

void setConstants(int dimx, int dimy, int dimz, int tdim) {
	// used in snow sim to convert snow particle
	// positions to wind grid space
	float temp[3] = {
	    (float)dimx / (float)SCENE_X,
	    (float)dimy / (float)SCENE_Y,
	    (float)dimz / (float)SCENE_Z
	};
	cudaMemcpyToSymbol(convert, temp, 3*sizeof(float), 
        0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(terrain_dim, &tdim, sizeof(int), 
        0, cudaMemcpyHostToDevice);
}

void initGPUWindResources(cudaArray **windVelArrayDevice, 
    int dimx, int dimy, int dimz, int tdim) {

    cudaChannelFormatDesc desc = cudaCreateChannelDesc<float4>();
    cudaExtent extent = make_cudaExtent(dimx, dimz, dimy);
    CUDA_SAFE_CALL(cudaMalloc3DArray(windVelArrayDevice, &desc, extent));

    wind_vel_tex.filterMode = cudaFilterModeLinear;
    wind_vel_tex.addressMode[0] = cudaAddressModeClamp;
    wind_vel_tex.addressMode[1] = cudaAddressModeClamp;
    wind_vel_tex.addressMode[2] = cudaAddressModeClamp;
    CUDA_SAFE_CALL(cudaBindTextureToArray(wind_vel_tex, *windVelArrayDevice));

	setConstants(dimx, dimy, dimz, tdim);
}

void windToGPU(float4 *windVel, cudaArray *windVelArrayDevice, 
    int dimx, int dimy, int dimz) {

    cudaMemcpy3DParms parm = {0};
    parm.srcPtr = make_cudaPitchedPtr(windVel, dimx*sizeof(float4), dimx, dimz);
    parm.dstArray = windVelArrayDevice;
    parm.extent = make_cudaExtent(dimx, dimz, dimy);
    parm.kind = cudaMemcpyHostToDevice;
    CUDA_SAFE_CALL(cudaMemcpy3D(&parm));
}

void freeGPUWindResources(cudaArray *windVelArrayDevice) {
    CUDA_SAFE_CALL(cudaUnbindTexture(wind_vel_tex));
    CUDA_SAFE_CALL(cudaFreeArray(windVelArrayDevice));
}

void getTerrainMapFromGPU(float4 *terrainMap, int tdim) {
	float4* device_map = get_terrain_vertices();
	CUDA_SAFE_CALL(cudaMemcpy(terrainMap, device_map, 
        sizeof(float4)*tdim*tdim, cudaMemcpyDeviceToHost));
}