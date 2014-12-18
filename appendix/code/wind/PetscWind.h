#ifndef PETSC_WIND_H
#define PETSC_WIND_H

#include <vector_types.h>
#include <vector>
#include <GL/glew.h>

#include "Camera.h"
#include "ShaderProgram.h"
#include "Wind.h"
#include "petscpp/petscppksp.h"

struct wind_shader_t {
	ShaderProgram pressure;
	ShaderProgram obstacle;
	ShaderProgram velocity;
};

class PetscWind : public Wind {
public:
	PetscWind();
	~PetscWind();
	
	void simulateWind() override;
	void updateObstacles() override;
	void renderObstacles(Camera *camera) override;
	void renderPressure(Camera *camera) override;
	void renderVelocityLines(Camera *camera) override;
	
private:
	// Dimensions of the wind velocity field, not including the boundary
	int3 dim;
	
	// Total number of unknowns, including the boundary
	int N;
	
	// Width and height of the terrain map
	int tdim;
	
	// Stores the matrix for the poisson problem
	petscpp::Matrix A;
	
	// Stores the pressure vector, result of solving the poisson problem
	petscpp::Vector p;
	
	// Stores the right hand side for the poisson problem
	petscpp::Vector b;
	
	// Solver for the poisson problem
	petscpp::KspSolver solver;
	
	// Vector storing the terrain map for obstacle creation
	std::vector<float4> map;
	
	// Velocity field
	std::vector<float4> velocity, prevVelocity, velocityStar;
	
	// Vector of bit-masks describing the occupation of obstacles 
	std::vector<int> obstacles;
	
	// cudaArray storing the wind velocity on the device
	cudaArray* windVelArrayDevice;
	
	// Struct storing the wind visualization shaders
	wind_shader_t shaders;
	
	// Vertex buffer object for visualization of wind simulator data
	GLuint point_vbo;
	
	// Self-advection, see Saltvik's or Eidissen's thesis
	void advect(double dt);
	
	// Projection, see Saltvik's or Eidissen's thesis
	void project(double dt);
	
	// Sets the velocities on the boundary of the domain
	void setVelocityBoundaries(std::vector<float4> &velocity);
	
	// Sets the pressure on the boundary of the domain
	void setBoundaries(petscpp::Vector &vec);
	
	// Calculate the global index of a 3d coordinate
	int getIndex(int x, int y, int z);
	
	// Helper function until a 3D volume class can be created
	template <class T>
	T get(const std::vector<T> &data, int x, int y, int z);
	
	// Helper function until a 3D volume class can be created
	template <class T>
	T get(T const * const data, int x, int y, int z);
	
	// Function for trilinear interpolation for the advection
	float4 trilinearInterpolation(const std::vector<float4> &data, 
		float x, float y, float z);
	
	void initVisualization();
	
	void setupSolution(double dt);
	
	void setupMatrix();
};

#endif
