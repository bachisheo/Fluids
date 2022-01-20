#pragma once
#define IX(x, y, z) ((x) + (y) * N + (z) * N * N)
#include "FluidLogic.h"
using std::vector;
class FluidCube {
public:
	FluidCube(int size, float diffusion, float viscosity, float dt, int iterCount);
	void addDensity(int x, int y, float value);
	void addDensity(vector<vector<float>>& x, const vector<vector<float>> & x0);
	void update() ;
	void dens_step();
	void vel_step();
	void reset() ;
	void render(int pixelOnFluidParticle, sf::Image& image) ;
	~FluidCube() ;
	void diffuseVelocity() ;
	void diffuseDensity() ;
	void addVelocity(int x0, int y0, int x1, int y1, float density_value);
	void diffuse(vector<vector<float>>& x, vector<vector<float>> const& x0, float diff, int b);
	void advection(vector<vector<float>>& d, vector<vector<float>>& d0, vector<vector<float>>& _u, vector<vector<float>>& _v, float dt, int b) ;
	void advection(vector<vector<float>>& d, vector<vector<float>>& d0, float dt, int b) ;


	void projectVelocity() ;
	void projectVelocity(vector<vector<float>>& div, vector<vector<float>>& p);
	virtual void setBoundary(vector<vector<float>> & x, int b) ;
private:
	float dt;
	float visc;
	float _diff;
	int height, width, N;
	int solverIterations = 50;
	vector<vector<float>> s;
	vector<vector<float>> dens;
	vector<vector<float>> dens_prev;

	vector<vector<float>> u;
	vector<vector<float>> v;
	vector<vector<float>> u_prev;
	vector<vector<float>> v_prev;
};
void Gauss(const vector<vector<float>>& x, const vector<vector<float>> & x0);

