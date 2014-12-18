#include "PetscWind.h"
#include "config.h"
#include "PetscWind.cuh"
#include "defines.h"
#include "ParticleSystem.h"

#include <algorithm>
#include <iostream>

#define VOXEL_MASK_SELF		(1)
#define VOXEL_MASK_RIGHT	(2)
#define VOXEL_MASK_LEFT 	(4)
#define VOXEL_MASK_BOTTOM	(8)
#define VOXEL_MASK_TOP		(16)
#define VOXEL_MASK_FRONT	(32)
#define VOXEL_MASK_BACK		(64)

/* 
 * To use petsc, the following commandline arguments must be specified:
 * -vec_type, -mat_type, -ksp_type, -pc_type
 *
 * Example of usage:
 * ./snow -vec_type seq -mat_type aij -ksp_type cg -pc_type jacobi
 *
 * To run on the gpu:
 * ./snow -vec_type seqcusp -mat_type aijcusp -ksp_type cg -pc_type jacobi
 */

PetscWind::PetscWind() : 
	N((conf.wind_x+2) * (conf.wind_y+2) * (conf.wind_z+2)),
	obstacles(N), velocity(N), prevVelocity(N), velocityStar(N),
	tdim(conf.map_size / conf.step_size), map(tdim*tdim),
	A(N, N), p(N), b(N)
{
	dim = make_int3(conf.wind_x, conf.wind_y, conf.wind_z);
	
	setupMatrix();
	p.fill(0);
	b.fill(0);
	
	A.assemblyBegin();
	p.assemblyBegin();
	b.assemblyBegin();
	A.assemblyEnd();
	p.assemblyEnd();
	b.assemblyEnd();
	
	A.setNullSpace();
	
	solver.setOperators(A, A);
	solver.createPreconditioner();
	solver.setTolerances(1E-3, 1E-3, 1E6, 15);
	solver.setup();
	
	for (int i = 0; i < N; i++) {
		velocity[i] = prevVelocity[i] = velocityStar[i] = {0};
		obstacles[i] = 0;
	}
	
	float4 v = make_float4(
		conf.wind_direction[0]*conf.wind_speed,
		conf.wind_direction[1]*conf.wind_speed,
		conf.wind_direction[2]*conf.wind_speed,
		0.0f
	);
	for (int z = 1; z < dim.z+1; z++) {
		for (int y = 1; y < dim.y+1; y++) {
			for (int x = 1; x < dim.x+1; x++) {
				int id = getIndex(x, y, z);
				velocity[id] = v;
			}
		}
	}
	
	initGPUWindResources(&windVelArrayDevice, 
		dim.x+2, dim.y+2, dim.z+2, tdim);
	windToGPU(prevVelocity.data(), windVelArrayDevice, 
		dim.x+2, dim.y+2, dim.z+2);
	
	initVisualization();
}

PetscWind::~PetscWind() {
	freeGPUWindResources(windVelArrayDevice);
}

void PetscWind::initVisualization() {
}

void PetscWind::renderObstacles(Camera *camera) {
}

void PetscWind::renderPressure(Camera *camera) {
}

void PetscWind::renderVelocityLines(Camera *camera) {
}

void PetscWind::simulateWind() {
	double const dt = conf.time_step;
	advect(dt);
	setupSolution(dt);
	p.fill(0);
	solver.solve(b, p);
	solver.printDebug();
	project(dt);
	windToGPU(prevVelocity.data(), windVelArrayDevice, 
		dim.x+2, dim.y+2, dim.z+2);
}

template <class T>
const T& clamp(const T& value, const T& min_value, const T& max_value) {
	return std::max(std::min(value, max_value), min_value);
}

float4 PetscWind::trilinearInterpolation(const std::vector<float4> &data, 
	float x, float y, float z) {

	float4 ret = {0};
	
	int i0 = (int) x;
	int i1 = i0 + 1;
	int j0 = (int) y;
	int j1 = j0 + 1;
	int k0 = (int) z;
	int k1 = k0 + 1;
	
	float s1 = x - i0;
	float s0 = 1 - s1;
	float t1 = y - j0;
	float t0 = 1 - t1;
	float u1 = z - k0;
	float u0 = 1 - u1;
	
	float lower = 0.0f, upper = 0.0f;
	
	// interpolate x
	lower = s0 * (t0 * data[getIndex(i0, j0, k0)].x 
		+ t1 * data[getIndex(i0, j1, k0)].x)
		+ s1 * (t0 * data[getIndex(i1, j0, k0)].x 
		+ t1 * data[getIndex(i1, j1, k0)].x);
	upper = s0 * (t0 * data[getIndex(i0, j0, k1)].x 
		+ t1 * data[getIndex(i0, j1, k1)].x)
		+ s1 * (t0 * data[getIndex(i1, j0, k1)].x 
		+ t1 * data[getIndex(i1, j1, k1)].x);
	ret.x = u0 * lower + u1 * upper;

	// interpolate y
	lower = s0 * (t0 * data[getIndex(i0, j0, k0)].y 
		+ t1 * data[getIndex(i0, j1, k0)].y)
		+ s1 * (t0 * data[getIndex(i1, j0, k0)].y 
		+ t1 * data[getIndex(i1, j1, k0)].y);
	upper = s0 * (t0 * data[getIndex(i0, j0, k1)].y 
		+ t1 * data[getIndex(i0, j1, k1)].y)
		+ s1 * (t0 * data[getIndex(i1, j0, k1)].y 
		+ t1 * data[getIndex(i1, j1, k1)].y);
	ret.y = u0 * lower + u1 * upper;

	// interpolate z
	lower = s0 * (t0 * data[getIndex(i0, j0, k0)].z 
		+ t1 * data[getIndex(i0, j1, k0)].z)
		+ s1 * (t0 * data[getIndex(i1, j0, k0)].z 
		+ t1 * data[getIndex(i1, j1, k0)].z);
	upper = s0 * (t0 * data[getIndex(i0, j0, k1)].z 
		+ t1 * data[getIndex(i0, j1, k1)].z)
		+ s1 * (t0 * data[getIndex(i1, j0, k1)].z 
		+ t1 * data[getIndex(i1, j1, k1)].z);
	ret.z = u0 * lower + u1 * upper;
	
	return ret;
}

void PetscWind::advect(double dt) {
	for (int z = 1; z < dim.z+1; z++) {
		for (int y = 1; y < dim.y+1; y++) {
			for (int x = 1; x < dim.x+1; x++) {
				int id = getIndex(x, y, z);
				
				float4 v = velocity[id];
				float xf = x - dt*v.x;
				float yf = y - dt*v.y;
				float zf = z - dt*v.z;
				
				// if outside of domain, clamp values
				xf = clamp(xf, 0.5f, dim.x+0.5f);
				yf = clamp(yf, 0.5f, dim.y+0.5f);
				zf = clamp(zf, 0.5f, dim.z+0.5f);
				
				velocityStar[id] = trilinearInterpolation(velocity, xf, yf, zf);
			}
		}
	}
	
	setVelocityBoundaries(velocityStar);
}

void PetscWind::project(double dt) {
	float h = 1.0f;
	float factor = 0.5f * dt / h;
	
	setBoundaries(p);
	
	PetscScalar * pressure = p.getArray();
	
	for (int z = 1; z < dim.z+1; z++) {
		for (int y = 1; y < dim.y+1; y++) {
			for (int x = 1; x < dim.x+1; x++) {
				int id = getIndex(x, y, z);
				
				prevVelocity[id] = velocity[id];
				
				float4 star = velocityStar[id];
				star.x -= factor * (get(pressure, x+1, y  , z  ) 
					- get(pressure, x-1, y  , z  ));
				star.y -= factor * (get(pressure, x  , y+1, z  ) 
					- get(pressure, x  , y-1, z  ));
				star.z -= factor * (get(pressure, x  , y  , z+1) 
					- get(pressure, x  , y  , z-1));
				
				velocity[id] = star;
			}
		}
	}
	
	p.restoreArray(&pressure);
	
	setVelocityBoundaries(velocity);
	setVelocityBoundaries(prevVelocity);
	
}

void PetscWind::setVelocityBoundaries(std::vector<float4> &velocity) {
	for (int z = 1; z < dim.z+1; z++) {
		for (int y = 1; y < dim.y+1; y++) {
			for (int x = 1; x < dim.x+1; x++) {
				
				const int id = getIndex(x, y, z);
				const int mask = obstacles[id];
				
				if (mask & VOXEL_MASK_SELF) {
					if (!(mask & VOXEL_MASK_LEFT)) {
						float4 old = get(velocity, x-1, y, z);
						velocity[id] = make_float4(0.0f, old.y, old.z, 0.0f);
						
					} else if (!(mask & VOXEL_MASK_RIGHT)) {
						float4 old = get(velocity, x+1, y, z);
						velocity[id] = make_float4(0.0f, old.y, old.z, 0.0f);
						
					} else if (!(mask & VOXEL_MASK_TOP)) {
						float4 old = get(velocity, x, y, z+1);
						velocity[id] = make_float4(old.x, old.y, 0.0f, 0.0f);
						
					} else if (!(mask & VOXEL_MASK_BOTTOM)) {
						float4 old = get(velocity, x, y, z-1);
						velocity[id] = make_float4(old.x, old.y, 0.0f, 0.0f);
						
					} else if (!(mask & VOXEL_MASK_FRONT)) {
						float4 old = get(velocity, x, y+1, z);
						velocity[id] = make_float4(old.x, 0.0f, old.z, 0.0f);
						
					} else if (!(mask & VOXEL_MASK_BACK)) {
						float4 old = get(velocity, x, y-1, z);
						velocity[id] = make_float4(old.x, 0.0f, old.z, 0.0f);
						
					} else {
						velocity[id] = {0.0f};
						
					}
				}
			}
		}
	}
	
	float4 value = {
		conf.wind_direction[0]*conf.wind_speed,
		conf.wind_direction[1]*conf.wind_speed,
		conf.wind_direction[2]*conf.wind_speed,
		0.0f
	};
	for (int z = 0; z < dim.z+2; z++) {
		for (int y = 0; y < dim.y+2; y++) {
			velocity[getIndex(0	  , y, z)] = value;
			velocity[getIndex(dim.x+1, y, z)] = value;
		}
	}
	
	for (int z = 0; z < dim.z+2; z++) {
		for (int x = 0; x < dim.x+2; x++) {
			velocity[getIndex(x, 0	  , z)] = value;
			velocity[getIndex(x, dim.y+1, z)] = value;
		}
	}
	
	for (int y = 0; y < dim.y+2; y++) {
		for (int x = 0; x < dim.x+2; x++) {
			velocity[getIndex(x, y, 0	  )] = value;
			velocity[getIndex(x, y, dim.z+1)] = value;
		}
	}
}

void PetscWind::setBoundaries(petscpp::Vector &vec) {
	PetscScalar * f = vec.getArray();
	for (int z = 0; z < dim.z+2; z++) {
		for (int y = 0; y < dim.y+2; y++) {
			f[getIndex(0	  , y, z)] = f[getIndex(1	, y, z)];
			f[getIndex(dim.x+1, y, z)] = f[getIndex(dim.x, y, z)];
		}
	}
	
	for (int z = 0; z < dim.z+2; z++) {
		for (int x = 0; x < dim.x+2; x++) {
			f[getIndex(x, 0	  , z)] = f[getIndex(x, 1	, z)];
			f[getIndex(x, dim.y+1, z)] = f[getIndex(x, dim.y, z)];
		}
	}
	
	for (int y = 0; y < dim.y; y++) {
		for (int x = 0; x < dim.x+2; x++) {
			f[getIndex(x, y, 0	  )] = f[getIndex(x, y, 1	)];
			f[getIndex(x, y, dim.z+1)] = f[getIndex(x, y, dim.z)];
		}
	}
	
	for (int z = 1; z < dim.z+1; z++) {
		for (int y = 1; y < dim.y+1; y++) {
			for (int x = 1; x < dim.x+1; x++) {
				
				const int id = getIndex(x, y, z);
				const int mask = obstacles[id];
				
				if (mask & VOXEL_MASK_SELF) {
					if (!(mask & VOXEL_MASK_LEFT)) {
						f[id] = get(f, x-1, y, z);
					} else if (!(mask & VOXEL_MASK_RIGHT)) {
						f[id] = get(f, x+1, y, z);
					} else if (!(mask & VOXEL_MASK_TOP)) {
						f[id] = get(f, x, y, z+1);
					} else if (!(mask & VOXEL_MASK_BOTTOM)) {
						f[id] = get(f, x, y, z-1);
					} else if (!(mask & VOXEL_MASK_FRONT)) {
						f[id] = get(f, x, y-1, z);
					} else if (!(mask & VOXEL_MASK_BACK)) {
						f[id] = get(f, x, y+1, z);
					} else {
						f[id] = 0.0f;
					}
				}
			}
		}
	}
	vec.restoreArray(&f);
}

void PetscWind::updateObstacles() {
	
	getTerrainMapFromGPU(map.data(), tdim);
	
	for (int z = 1; z < dim.z+1; z++) {
		for (int y = 1; y < dim.y+1; y++) {
			for (int x = 1; x < dim.x+1; x++) {
				float h = SCENE_Y * static_cast<float>(y) / dim.y;
				int sx = ((x-1) * tdim) / dim.x;
				int sz = ((z-1) * tdim) / dim.z;
				
				float4 const &v = map[(tdim*sz)+sx];
				obstacles[getIndex(x, y, z)] = ((h < v.y + v.w) & 1) 
					* VOXEL_MASK_SELF;
			}
		}
	}
	
	for (int z = 1; z < dim.z+1; z++) {
		for (int y = 1; y < dim.y+1; y++) {
			for (int x = 1; x < dim.x+1; x++) {
				obstacles[getIndex(x,y,z)] += 
						(get(obstacles, x-1, y, z) & 1) * VOXEL_MASK_LEFT +
						(get(obstacles, x+1, y, z) & 1) * VOXEL_MASK_RIGHT +
						(get(obstacles, x, y-1, z) & 1) * VOXEL_MASK_FRONT +
						(get(obstacles, x, y+1, z) & 1) * VOXEL_MASK_BACK +
						(get(obstacles, x, y, z-1) & 1) * VOXEL_MASK_BOTTOM +
						(get(obstacles, x, y, z+1) & 1) * VOXEL_MASK_TOP;
			}
		}
	}
}

int PetscWind::getIndex(int x, int y, int z) {
	return x + z*(dim.x+2) + y*(dim.x+2)*(dim.z+2);
}

template <class T>
T PetscWind::get(const std::vector<T> &data, int x, int y, int z) {
	return data[getIndex(x, y, z)];
}

template <class T>
T PetscWind::get(T const * const data, int x, int y, int z) {
	return data[getIndex(x, y, z)];
}

void PetscWind::setupSolution(double dt) {
	PetscScalar hx = 1 / static_cast<PetscScalar>(dim.x+1);
	PetscScalar hy = 1 / static_cast<PetscScalar>(dim.y+1);
	PetscScalar hz = 1 / static_cast<PetscScalar>(dim.z+1);
	
	PetscScalar hyz = hy*hy*hz*hz;
	PetscScalar hxz = hx*hx*hz*hz;
	PetscScalar hxy = hx*hx*hy*hy;
	
	PetscScalar * solution = b.getArray();
	
	for (int z = 1; z < dim.z+1; z++) {
		for (int y = 1; y < dim.y+1; y++) {
			for (int x = 1; x < dim.x+1; x++) {
				int id = getIndex(x, y, z);
				
				PetscScalar ux = (velocityStar[getIndex(x+1, y, z)].x 
					- velocityStar[getIndex(x, y, z)].x);
				PetscScalar uy = (velocityStar[getIndex(x, y+1, z)].y 
					- velocityStar[getIndex(x, y, z)].y);
				PetscScalar uz = (velocityStar[getIndex(x, y, z+1)].z 
					- velocityStar[getIndex(x, y, z)].z);
				
				solution[id] = ((obstacles[id] ^ 1) & 1) * (hx*hyz*ux 
															+ hy*hxz*uy 
															+ hz*hxy*uz);
			}
		}
	}
	
	b.restoreArray(&solution);
	setBoundaries(b);
}

void PetscWind::setupMatrix() {
	PetscScalar hx = 1 / static_cast<PetscScalar>(dim.x+1);
	PetscScalar hy = 1 / static_cast<PetscScalar>(dim.y+1);
	PetscScalar hz = 1 / static_cast<PetscScalar>(dim.z+1);
	
	PetscScalar hyz = hy*hy*hz*hz;
	PetscScalar hxz = hx*hx*hz*hz;
	PetscScalar hxy = hx*hx*hy*hy;
	
	for (int z = 1; z < dim.z+1; z++) {
		for (int y = 1; y < dim.y+1; y++) {
			for (int x = 1; x < dim.x+1; x++) {
				PetscInt row = getIndex(x, y, z);
				A.set(row, getIndex(x-1, y, z), hyz);
				A.set(row, getIndex(x+1, y, z), hyz);
				A.set(row, getIndex(x, y-1, z), hxz);
				A.set(row, getIndex(x, y+1, z), hxz);
				A.set(row, getIndex(x, y, z-1), hxy);
				A.set(row, getIndex(x, y, z+1), hxy);
				A.set(row, row, -2*(hyz+hxz+hxy));
			}
		}
	}
}
