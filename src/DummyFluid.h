#pragma once
#define IX(x, y, z) ((x) + (y) * N + (z) * N * N)
#include "FluidLogic.h"
using std::vector;
class FluidCube {
public:
	FluidCube(int size, float diffusion, float viscosity, int iterCount);
	void addDensity(int x, int y, float value);
	void addDensity(vector<vector<float>>& x, const vector<vector<float>>& x0, float dt);
	void addVelocity(int x, int y, int dx, int dy, float dt);
	void update(float dt);
	void dens_step(float dt);
	void vel_step(float dt);
	void reset();
	void render(int pixelOnFluidParticle, sf::Image& image);
	void diffuse(vector<vector<float>>& x, vector<vector<float>> const& x0, float diff, int b, float dt);
	void advection(vector<vector<float>>& d, vector<vector<float>>& d0, vector<vector<float>>& _u, vector<vector<float>>& _v, float dt, int b);
	void advection(vector<vector<float>>& d, vector<vector<float>>& d0, float dt, int b);

	float GetDiffusionCoef() const { return _diff; }
	float GetViscosityCoef() const { return visc; }
	int GetCursorRadiusCoef() const { return _cursorRadius; }
	void SetDiffusionCoef(float diff) { _diff = diff; }
	void SetViscosityCoef(float visc) { this->visc = visc; }
	void SetCursorRadiusCoef(int cursorRadius) { _cursorRadius = cursorRadius; }
	float speed = 50;

	void projectVelocity(vector<vector<float>>& div, vector<vector<float>>& p);
	virtual void setBoundary(vector<vector<float>>& x, int b);
private:
	float visc;
	int _cursorRadius = 1;
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
void Gauss(const vector<vector<float>>& x, const vector<vector<float>>& x0, float koef = 1);

