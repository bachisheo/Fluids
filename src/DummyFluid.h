#pragma once
#define IX(x, y, z) ((x) + (y) * N + (z) * N * N)
#include "FluidLogic.h"

class FluidCube : public FluidLogic{
public:
	FluidCube(int size, int diffusion, int viscosity, float dt);
	void AddSourceOfDensity(int x0, int y0, int x, int y) override;
	void update() override;
	void reset() override;
	void render(int size, sf::Image& image) override;
	~FluidCube() override;
protected:
	void diffuseVelocity() override;
	void diffuseDensity() override;
	void advect() override;
	void projectVelocity() override;
	virtual void setBoundary(int d, Field x) override;
private:
	void AddDensity(int x, int y, int z, float amount);
	int size;
	float dt;
	float diff;
	float visc;

	float* s;
	float* density;

	float* Vx;
	float* Vy;
	float* Vz;

	float* Vx0;
	float* Vy0;
	float* Vz0;
};
