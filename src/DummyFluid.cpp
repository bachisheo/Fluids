#include "DummyFluid.h"

#include <corecrt_math.h>
#include <cstdlib>

FluidCube::FluidCube(int size, float diffusion, float viscosity, float dt, int iterCount) : solverIterations(iterCount)
{
	N = width = height = size;
	_diff = diffusion;
	visc = viscosity;
	this->dt = dt;
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
	dens_prev[x][y] = value;
}
void FluidCube::addDensity(vector<vector<float>>& x, const vector<vector<float>>& x0)
{
	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			x[i][j] += dt * x0[i][j];
}


FluidCube::~FluidCube()
{
}
void Gauss(vector<vector<float>>& x, const vector<vector<float>> & x0, float koef)
{
	for (int i = 1; i < x.size() - 1; i++) {
		for (int j = 1; j < x[i].size() - 1; j++) {
			x[i][j] = (x0[i][j] + koef * (x[i - 1][j] + x[i + 1][j] +
				x[i][j - 1] + x[i][j + 1])) / (4 * koef);
		}
	}
}

void FluidCube::addVelocity(int x0, int y0, int x1, int y1, float density_value) {
	int dx = x1 - x0;
	int dy = y1 - y0;
	if (abs(dx) >= abs(dy)) {
		int yi = 1;
		if (dy < 0) {
			yi = -1;
			dy = -dy;
		}
		int D = 2 * dy - dx;
		int y = y0;
		int elementID = 0;
		if (dx > 0) {
			for (int x = x0; x <= x1; x++) {
				elementID = x + (width + 2) * y;
				dens_prev[x][y] = 1;
				u_prev[x][y] = (x1 - x0) / dt;
				v_prev[x][y] = (y1 - y0) / dt;
				if (D > 0) {
					y += yi;
					D -= 2 * dx;
				}
				D += 2 * dy;
			}
		}
		else if (dx < 0) {
			D = 2 * dy + dx;
			for (int x = x0; x >= x1; x--) {
				//elementID = x + (xSize + 2) * y;
				dens_prev[x][y] = 1;
				u_prev[x][y] = (x1 - x0) / dt;
				v_prev[x][y] = (y1 - y0) / dt;
				if (D > 0) {
					y += yi;
					D += 2 * dx;
				}
				D += 2 * dy;
			}
		}
	}
	else if (abs(dy) > abs(dx)) {
		int xi = 1;
		if (dx < 0) {
			xi = -1;
			dx = -dx;
		}
		int D = 2 * dx - dy;
		int x = x0;
		int elementID = 0;
		if (dy > 0) {
			for (int y = y0; y <= y1; y++) {
				//elementID = x + (xSize + 2) * y;
				dens_prev[x][y] = 1;
				u_prev[x][y] = (x1 - x0) / dt;
				v_prev[x][y] = (y1 - y0) / dt;
				if (D > 0) {
					x += xi;
					D -= 2 * dy;
				}
				D += 2 * dx;
			}
		}
		else if (dy < 0) {
			D = 2 * dx + dy;
			for (int y = y0; y >= y1; y--) {
				//elementID = x + (xSize + 2) * y;
				dens_prev[x][y] = 1;
				u_prev[x][y] = (x1 - x0) / dt;
				v_prev[x][y] = (y1 - y0) / dt;
				if (D > 0) {
					x += xi;
					D += 2 * dy;
				}
				D += 2 * dx;
			}
		}
	}
}
//2. диффузия: плотность и скорость между соседними клетками меняется

void FluidCube::diffuse(vector<vector<float>>& x, vector<vector<float>> const& x0, float diff, int b)
{
	float a = dt * diff * width * height;
	for (int k = 0; k < solverIterations; k++) {
		for (int i = 1; i < x.size() - 1; i++) {
			for (int j = 1; j < x[i].size() - 1; j++) {
				x[i][j] = (x0[i][j] + a * (x[i - 1][j] + x[i + 1][j] +
					x[i][j - 1] + x[i][j + 1])) / ( 1 + 4 * a);
			}
		}
		setBoundary(x, b);
	}
}

void FluidCube::setBoundary(vector<vector<float>> & x, int b)
{
	float neg = -1.0;
	for (int i = 1; i < N; i++) {
		x[0][i] = x[1][i];
		x[N - 1][i] = x[N - 2][i];
		x[i][0] = x[i][1];
		x[i][N - 1] =  x[i][N - 2];


		/*if (b == 1)
		{
			//x[0][i]
		}
		
		x[0][i] = b == 1 ? neg * x[1][i] : x[1][i];
		x[N+1][i] = b == 1 ? neg * x[N][i] : x[N][i];
		x[i][0] = b == 2 ?  neg * x[i][1] : x[i][1];
		x[i][N+1] = b == 2 ? neg* x[i][N] : x[i][N];
		*/
	}
	x[0][0] = 0.5f * (x[1][0] + x[0][1]);
	x[0][N-1] = 0.5f * (x[1][N -2] + x[0][N - 1]);
	x[N-1][0] = 0.5f * (x[N-2][0] + x[N-1][1]);
	x[N-1][N-1] = 0.5f * (x[N-1][N-2] + x[N-2][N-1]);
}

//3. адвекция -- перенос потоков

void FluidCube::advection(vector<vector<float>> &d, vector<vector<float>> &d0, float dt, int b)
{
	advection(d, d0, u, v, dt, b);
}


void FluidCube::advection(vector<vector<float>> & d, vector<vector<float>>& d0, vector<vector<float>> &_u, vector<vector<float>>& _v, float dt, int b)
{
	float dt0 = dt * N;
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			float x = i - dt0 * _u[i][j];
			float y = j - dt0 * _v[i][j];
			if (x < 0.5)
				x = 0.5;
			if (x > N + 0.5)
				x = N + 0.5;
			//int i0 = std::min((int)x, height - 2);
			auto i0 = (int)x;
			int i1 = std::min(i0 + 1, height- 1);
			if (y < 0.5)
				y = 0.5;
			if (y > N + 0.5)
				y = N + 0.5;
			//int j0 = std::min((int)y, width - 2);
			auto j0 = (int)y;

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


void FluidCube::update()
{
	vel_step();
	dens_step();
	dens_prev = vector(height, vector<float>(width));
	u_prev = vector(height, vector<float>(width));
	v_prev = vector(height, vector<float>(width));
	//и рисуем денс
}

//шаги 1-3 можно объединить
//это решение плотности
void FluidCube::dens_step()
{
	addDensity(dens, dens_prev);
	std::swap(dens_prev, dens);
	diffuse( dens, dens_prev, _diff, 0);
	std::swap(dens_prev, dens);
	advection(dens, dens_prev, dt, 0);

}
void FluidCube::vel_step()
{
	addDensity(u, u_prev); 
	addDensity(v, v_prev);

	std::swap(u_prev, u);
	diffuse(u, u_prev, visc, 1);
	std::swap(v_prev, v);
	diffuse(v, v_prev, visc, 2);

	projectVelocity( u_prev, v_prev);

	std::swap(u_prev, u); 
	std::swap(v_prev, v);
	advection(u, u_prev, u_prev, v_prev, dt, 1);
	advection(v, v_prev, u_prev, v_prev, dt, 2);
	projectVelocity(u_prev, v_prev);

}
void FluidCube::projectVelocity(vector<vector<float>> & div, vector<vector<float>>& p)
{
	int i, j, k;
	float h;
	h = 1.0 / N;
	for (i = 1; i < N - 1; i++) {
		for (j = 1; j < N - 1; j++) {
			div[i][j] = -0.5 * h * (u[i + 1][j] - u[i-1][j] +
				v[i][j+1] - v[i][j-1]);
			p[i][j] = 0;
		}
	}
	setBoundary(div, 0);
	setBoundary(p, 0);
	for (k = 0; k < solverIterations; k++) {
		Gauss(p, div, 1);
		setBoundary(p, 0);
	}
	for (i = 1; i < height - 1; i++) {
		for (j = 1; j < width - 1; j++) {
			u[i][j] -= 0.5 * (p[i+1][j] - p[i-1][j]) / h;
			v[i][j] -= 0.5 * (p[i][j+1] - p[i][j-1]) / h;
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
	for (int i = 0; i < height ; i++) {
		for (int j = 0; j < width; j++) {
			if (dens[i][j] > 1) {
				dens[i][j] = 1;
			}
			for (int ii = 0; ii <= pixelOnFluidParticle; ii++) {
				for (int jj = 0; jj <= pixelOnFluidParticle; jj++) {
					auto a = 255 * dens[i][j];
					image.setPixel(i * pixelOnFluidParticle + ii, j * pixelOnFluidParticle + jj, sf::Color(a, a, a));
				}
			}
			//density0.values[i + j * (xSize + 2)] = density.values[i + j * (xSize + 2)];
			//u0.values[i + j * (xSize + 2)] = u.values[i + j * (xSize + 2)];
			//v0.values[i + j * (xSize + 2)] = v.values[i + j * (xSize + 2)];
		}
	}
}


