#include "DummyFluid.h"

#include <corecrt_math.h>
#include <cstdlib>
#include <iostream>

FluidCube::FluidCube(int size, float diffusion, float viscosity, int iterCount) : solverIterations(iterCount)
{
	N = width = height = size;
	_diff = diffusion;
	visc = viscosity;
	s = vector(height, vector<float>(width));
	dens = vector(height, vector<float>(width));
	dens_prev = vector(height, vector<float>(width));
	u = vector(height, vector<float>(width));
	v = vector(height, vector<float>(width));
	u_prev = vector(height, vector<float>(width));
	v_prev = vector(height, vector<float>(width));
}
///1. Инициализация
void FluidCube::addDensity(int x, int y, float value)
{
	x = std::clamp(x, 0, N - 1);
	y = std::clamp(y, 0, N - 1);
	int r = 5;
	for (int i = -r; i <= r; i++)
		for (int j = -r; j <= r; j++)
		{
			int x1 = x + i, y1 = y + j;
			if (x1 >= 0 && x1 < N && y1>0 && y1 < N)
				dens_prev[x1][y1] = value;
		}
}
void FluidCube::addDensity(vector<vector<float>>& x, const vector<vector<float>>& x0, float dt)
{
	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			x[i][j] += dt * x0[i][j];
}

void FluidCube::addVelocity(int x, int y, int dx, int dy, float dt)
{
	x = std::clamp(x, 0, height - 1);
	y = std::clamp(y, 0, width - 1);
	v_prev[x][y] = (float)dy / dt * speed;
	u_prev[x][y] = (float)dx / dt * speed;


}

void Gauss(vector<vector<float>>& x, const vector<vector<float>>& x0, float koef)
{
	for (int i = 1; i < x.size() - 1; i++) {
		for (int j = 1; j < x[i].size() - 1; j++) {
			x[i][j] = (x0[i][j] + koef * (x[i - 1][j] + x[i + 1][j] +
				x[i][j - 1] + x[i][j + 1])) / (4 * koef);
		}
	}
}


//2. диффузия: плотность и скорость между соседними клетками меняется

void FluidCube::diffuse(vector<vector<float>>& x, vector<vector<float>> const& x0, float diff, int b, float dt)
{
	float a = dt * diff;
	for (int k = 0; k < solverIterations; k++) {
		for (int i = 1; i < x.size() - 1; i++) {
			for (int j = 1; j < x[i].size() - 1; j++) {
				x[i][j] = (x0[i][j] + a * (x[i - 1][j] + x[i + 1][j] +
					x[i][j - 1] + x[i][j + 1])) / (1 + 4 * a);
			}
		}
		setBoundary(x, b);
	}
}

void FluidCube::setBoundary(vector<vector<float>>& x, int b)
{
	float neg = -1.0;
	for (int i = 1; i < N; i++) {
		x[0][i] = x[1][i];
		x[N - 1][i] = x[N - 2][i];
		x[i][0] = x[i][1];
		x[i][N - 1] = x[i][N - 2];


		if (b == 1)
		{
			x[0][i] *= neg;
			x[N - 1][i] *= neg;
		}
		if (b == 2)
		{
			x[i][0] *= neg;
			x[i][N - 1] *= neg;
		}

	}
	x[0][0] = 0.5f * (x[1][0] + x[0][1]);
	x[0][N - 1] = 0.5f * (x[1][N - 2] + x[0][N - 1]);
	x[N - 1][0] = 0.5f * (x[N - 2][0] + x[N - 1][1]);
	x[N - 1][N - 1] = 0.5f * (x[N - 1][N - 2] + x[N - 2][N - 1]);
}

//3. адвекция -- перенос потоков

void FluidCube::advection(vector<vector<float>>& d, vector<vector<float>>& d0, float dt, int b)
{
	advection(d, d0, u, v, dt, b);
}


void FluidCube::advection(vector<vector<float>>& d, vector<vector<float>>& d0, vector<vector<float>>& _u, vector<vector<float>>& _v, float dt, int b)
{
	float dt0 = dt;
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			float x = i - dt0 * _u[i][j];
			float y = j - dt0 * _v[i][j];

			if (x < 0.5)
				x = 0.5;
			if (x > N + 0.5)
				x = N + 0.5;
			int i0 = std::min((int)x, height - 2);

			int i1 = std::min(i0 + 1, height - 1);
			if (y < 0.5)
				y = 0.5;
			if (y > N + 0.5)
				y = N + 0.5;
			int j0 = std::min((int)y, width - 1);

			int j1 = std::min(j0 + 1, width - 1);
			float s1 = x - i0;
			float s0 = 1 - s1;
			float t1 = y - j0;
			float t0 = 1 - t1;
			d[i][j] = s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1]) +
				s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);
		}
	}
	setBoundary(d, b);
}


void FluidCube::update(float dt)
{
	
	vel_step(dt);
	dens_step(dt);
	dens_prev = vector(height, vector<float>(width));
	u_prev = vector(height, vector<float>(width));
	v_prev = vector(height, vector<float>(width));
	//и рисуем денс
}

//шаги 1-3 можно объединить
//это решение плотности
void FluidCube::dens_step(float dt)
{
	addDensity(dens, dens_prev, dt);
	std::swap(dens_prev, dens);
	diffuse(dens, dens_prev, _diff, 0, dt);
	std::swap(dens_prev, dens);
	advection(dens, dens_prev, dt, 0);

}
void FluidCube::vel_step(float dt)
{
	addDensity(u, u_prev, dt);
	addDensity(v, v_prev, dt);

	std::swap(u_prev, u);
	diffuse(u, u_prev, visc, 1, dt);
	std::swap(v_prev, v);
	diffuse(v, v_prev, visc, 2, dt);

	projectVelocity(u_prev, v_prev);

	std::swap(u_prev, u);
	std::swap(v_prev, v);

	advection(u, u_prev, u_prev, v_prev, dt, 1);
	advection(v, v_prev, u_prev, v_prev, dt, 2);
	projectVelocity(u_prev, v_prev);

}
void FluidCube::projectVelocity(vector<vector<float>>& p, vector<vector<float>>& div)
{
	int i, j, k;
	float h;
	h = 1.0 / N;
	for (i = 1; i < N - 1; i++) {
		for (j = 1; j < N - 1; j++) {
			div[i][j] = -0.5 * h * (u[i + 1][j] - u[i - 1][j] +
				v[i][j + 1] - v[i][j - 1]);
			p[i][j] = 0;
		}
	}
	setBoundary(div, 0);
	setBoundary(p, 0);
	for (k = 0; k < solverIterations; k++) {
		for (int i = 1; i < p.size() - 1; i++) {
			for (int j = 1; j < p[i].size() - 1; j++) {
				p[i][j] = (div[i][j] + p[i - 1][j] + p[i + 1][j] +
					p[i][j - 1] + p[i][j + 1]) / 4;
			}
		}
		setBoundary(p, 0);
	}
	for (i = 1; i < height - 1; i++) {
		for (j = 1; j < width - 1; j++) {
			u[i][j] -= 0.5 * (p[i + 1][j] - p[i - 1][j]) / h;
			v[i][j] -= 0.5 * (p[i][j + 1] - p[i][j - 1]) / h;
		}
	}
	setBoundary(u, 1);
	setBoundary(v, 2);
}

void FluidCube::reset()
{
	s = vector(height, vector<float>(width, 0));
	dens = vector(height, vector<float>(width, 0));
	dens_prev = vector(height, vector<float>(width, 0));
	u = vector(height, vector<float>(width, 0));
	v = vector(height, vector<float>(width, 0));
	u_prev = vector(height, vector<float>(width, 0));
	v_prev = vector(height, vector<float>(width, 0));
}

void FluidCube::render(int pixelOnFluidParticle, sf::Image& image)
{
	//draw pixel by pixel
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			auto coef = std::min(1.f, dens[i][j]);

			for (int ii = 0; ii <= pixelOnFluidParticle; ii++) {
				for (int jj = 0; jj <= pixelOnFluidParticle; jj++) {
					auto a = 255 * std::pow(coef,0.5f);

					image.setPixel(i * pixelOnFluidParticle + ii, j * pixelOnFluidParticle + jj, sf::Color(255*coef, 0, 255 * (1-coef), a));
				}
			}
			//density0.values[i + j * (xSize + 2)] = density.values[i + j * (xSize + 2)];
			//u0.values[i + j * (xSize + 2)] = u.values[i + j * (xSize + 2)];
			//v0.values[i + j * (xSize + 2)] = v.values[i + j * (xSize + 2)];
		}
	}
}


