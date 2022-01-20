#pragma once
#include "FluidLogic.h"
#include <SFML/Graphics.hpp>

class Fluid : FluidLogic
{
public:
	//constructor and destructor
	Fluid();
	Fluid(int f_xSize, int f_ySize, float f_diff, float f_visc, float f_dt, int f_solverIterations);
	~Fluid() override;

	//the public functions
	void AddSourceOfDensity(int x0, int y0, int x, int y) override;
	void update() override;
	void reset() override;
	void render(int size, sf::Image &image) override;

private:
	//Field property declarations
	Field u;
	Field v;
	Field density;
	Field u0;
	Field v0;
	Field density0;
	Field uSource;
	Field vSource;
	Field densitySource;
	Field div;
	Field p;

	//Physical Property declarations
	float diff;
	float visc;

	//Simulation parameters
	float dt;
	int solverIterations;
	int xSize;
	int ySize;

	//fluid physics functions
	void diffuseVelocity() override;
	void diffuseDensity() override;
	void advection() override;
	void advectVelocity() ;
	void advectDensity() ;
	void projectVelocity() override;
	void setBoundary(int d, Field x) override;
	void swapPointers(float *&x0, float *&x1);
};

