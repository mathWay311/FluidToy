
#include <iostream>
#include "fluid.h"
#include <cmath>


using enum debug_mode;

FluidDomain::FluidDomain(unsigned int numx, unsigned numy, debug_mode mode_)
{
	enum field_type {U_FIELD, V_FIELD, S_FIELD};
	double overRelaxation = 1.9;
	density = 1000;
	mode = mode_;
	iteration_counter = 0;
	num_cells_x = numx + 2;
	num_cells_y = numy + 2;
	num_cells = num_cells_y * num_cells_x;
	U.resize(num_cells);
	V.resize(num_cells);
	newU.resize(num_cells);
	newV.resize(num_cells);
	P.resize(num_cells);
	S.resize(num_cells);
	M.resize(num_cells);
	newM.resize(num_cells);

	if (mode == FULL)
	{
		std::cout << "Total number of cells: " << U.size() << std::endl;
	}
}

FluidDomain::FluidDomain()
{
	enum field_type {U_FIELD, V_FIELD, S_FIELD};
	double overRelaxation = 1.9;
	density = 1000;
	mode = FULL;
	iteration_counter = 0;
	num_cells_x = 0;
	num_cells_y = 0;
	num_cells = num_cells_y * num_cells_x;

	U.resize(num_cells);
	V.resize(num_cells);
	newU.resize(num_cells);
	newV.resize(num_cells);
	P.resize(num_cells);
	S.resize(num_cells);
	M.resize(num_cells);
	newM.resize(num_cells);

	if (mode == FULL)
	{
		std::cout << "Total number of cells: " << U.size() << std::endl;
	}
}

void FluidDomain::integrate(double dt, double gravity)
{
	int n = num_cells_y;
	for (unsigned int i = 1; i < num_cells_x; i++) 
	{
		for (unsigned int j = 1; j < num_cells_y - 1; j++) 
		{
			if (S[i*n + j] != 0.0 && S[i*n + j-1] != 0.0) V[i*n + j] += gravity * dt;
		}	 
	}
}


void FluidDomain::solveIncompressibility(unsigned int numIters, double dt) 
{
			unsigned int n = num_cells_y;
			double cp = density * num_cells_y / dt;

			for (unsigned int iter {0}; iter < numIters; iter++) 
			{
				for (unsigned int i {1}; i < num_cells_x - 1; i++) 
				{
					for (unsigned int j {1}; j < num_cells_y - 1; j++) 
					{

						if (S[i*n + j] == 0.0)
							continue;

						double sx0 = S[(i-1)*n + j];
						double sx1 = S[(i+1)*n + j];
						double sy0 = S[i*n + j-1];
						double sy1 = S[i*n + j+1];
						double s = sx0 + sx1 + sy0 + sy1;
						if (s == 0.0)
							continue;

						double div = U[(i+1)*n + j] - U[i*n + j] + V[i*n + j+1] - V[i*n + j];

						double p = -div / s;
						p *= overRelaxation;
						P[i*n + j] += cp * p;

						U[i*n + j] -= sx0 * p;
						U[(i+1)*n + j] += sx1 * p;
						V[i*n + j] -= sy0 * p;
						V[i*n + j+1] += sy1 * p;
					}
				}
			}
		}

void FluidDomain::extrapolate() {
			int n = num_cells_y;
			for (unsigned int i = 0; i < num_cells_x; i++) {
				U[i*n + 0] = U[i*n + 1];
				U[i*n + num_cells_y - 1] = U[i*n + num_cells_y - 2]; 
			}
			for (unsigned int j = 0; j < num_cells_y; j++) {
				V[0*n + j] = V[1*n + j];
				V[(num_cells_x - 1)*n + j] = V[(num_cells_x - 2)*n + j];
			}
		}

double FluidDomain::sampleField(double x, double y, field_type field) {
			double n = num_cells_y;
			double h = num_cells_y;
			double h1 = 1.0 / h;
			double h2 = 0.5 * h;

			x = std::max(std::min(x, num_cells_x * h), h);
			y = std::max(std::min(y, num_cells_y * h), h);

			double dx = 0.0;
			double dy = 0.0;

			switch (field) {
				case U_FIELD: dy = h2; break;
				case V_FIELD: dx = h2; break;
				case S_FIELD: dx = h2; dy = h2; break;
				default:break;

			}

			double x0 = std::min(std::floor((x-dx)*h1), (double)num_cells_x - 1);
			double tx = ((x-dx) - x0*h) * h1;
			double x1 = std::min(x0 + 1, (double)num_cells_x - 1);
			
			double y0 = std::min(std::floor((y-dy)*h1), (double)num_cells_y - 1);
			double ty = ((y-dy) - y0*h) * h1;
			double y1 = std::min(y0 + 1, (double)num_cells_y - 1);

			double sx = 1.0 - tx;
			double sy = 1.0 - ty;

			double val;
			switch (field) {
				case U_FIELD:
					val = sx*sy * U[x0*n + y0] +
					tx*sy * U[x1*n + y0] +
					tx*ty * U[x1*n + y1] +
					sx*ty * U[x0*n + y1];
					break;
				case V_FIELD: 
					val = sx*sy * V[x0*n + y0] +
					tx*sy * V[x1*n + y0] +
					tx*ty * V[x1*n + y1] +
					sx*ty * V[x0*n + y1];
					break;
				case S_FIELD: 
					val = sx*sy * M[x0*n + y0] +
					tx*sy * M[x1*n + y0] +
					tx*ty * M[x1*n + y1] +
					sx*ty * M[x0*n + y1];
					break;
				default:
					break;
			}

			
			return val;
		}

double FluidDomain::avgU(unsigned int i, unsigned int j) 
{
	int n = num_cells_y;
	U[i*n + j] = (U[i*n + j-1] + U[i*n + j] +
		U[(i+1)*n + j-1] + U[(i+1)*n + j]) * 0.25;
	return U[i*n + j];
}

double FluidDomain::avgV(unsigned int i, unsigned int j) 
{
	int n = num_cells_y;
	V[i*n + j] = (V[i*n + j-1] + V[i*n + j] +
		V[(i+1)*n + j-1] + V[(i+1)*n + j]) * 0.25;
	return V[i*n + j];
}

void FluidDomain::advectVel(double dt) {

			newU = U;
			newV = V;

			unsigned int n = num_cells_y;
			double h = num_cells_y;
			double h2 = 0.5 * h;

			for (unsigned int i = 1; i < num_cells_x; i++) {
				for (unsigned int j = 1; j < num_cells_y; j++) {

					//cnt++;

					// u component
					if (S[i*n + j] != 0.0 && S[(i-1)*n + j] != 0.0 && j < num_cells_y - 1) {
						double x = i*h;
						double y = j*h + h2;
						double u = U[i*n + j];
						double v = avgV(i, j);
//						var v = sampleField(x,y, V_FIELD);
						x = x - dt*u;
						y = y - dt*v;
						u = sampleField(x,y, U_FIELD);
						newU[i*n + j] = u;
					}
					// v component
					if (S[i*n + j] != 0.0 && S[i*n + j-1] != 0.0 && i < num_cells_x - 1) {
						double x = i*h + h2;
						int y = j*h;
						double u = avgU(i, j);
//						var u = this.sampleField(x,y, U_FIELD);
						double v = V[i*n + j];
						x = x - dt*u;
						y = y - dt*v;
						v = sampleField(x,y, V_FIELD);
						newV[i*n + j] = v;
					}
				}	 
			}

			U = newU;
			V = newV;
		}

void FluidDomain::advectSmoke(double dt) 
{

			newM = M;

			double n = num_cells_y;
			double h = num_cells_y;
			double h2 = 0.5 * h;

			for (unsigned int i = 1; i < num_cells_x - 1; i++) {
				for (unsigned int j = 1; j < num_cells_y - 1; j++) {

					if (S[i*n + j] != 0.0) {
						double u = (U[i*n + j] + U[(i+1)*n + j]) * 0.5;
						double v = (V[i*n + j] + V[i*n + j+1]) * 0.5;
						double x = i*h + h2 - dt*u;
						double y = j*h + h2 - dt*v;

						newM[i*n + j] = sampleField(x,y, S_FIELD);
 					}
				}	 
			}
			M = newM;
}

void FluidDomain::simulate(double dt, double gravity, unsigned int numIters) 
{
	if (mode == MIN || mode == FULL) print_msg("=== Начало итерации " + std::to_string(iteration_counter) + " ===");
	if (mode == MIN || mode == FULL) print_msg("Начало симуляции...");
	if (mode == MIN || mode == FULL) print_msg("Интегрирование...");
	integrate(dt, gravity);
	if (mode == MIN || mode == FULL) print_msg("Решение несжимаемости...");
	solveIncompressibility(numIters, dt);
	//if (mode == MIN || mode == FULL) print_msg("Экстраполяция...");
	extrapolate();
	if (mode == MIN || mode == FULL) print_msg("Адвекция скорости...");
	advectVel(dt);
	if (mode == MIN || mode == FULL) print_msg("Адвекция дыма...");
	advectSmoke(dt);
	if (mode == MIN || mode == FULL) print_msg("=== Конец итерации " + std::to_string(iteration_counter) + " ===");
	iteration_counter++;

}


void FluidDomain::debug_print_field(field_type type)
{
	for (unsigned int i = 0; i < num_cells_x; i++) 
	{
		for (unsigned int j = 0; j < num_cells_y; j++)
		{
			switch(type)
			{
				case U_FIELD:
					std::cout << U[i*num_cells_y + j] << " ";
					break;
				case V_FIELD:
					std::cout << V[i*num_cells_y + j] << " ";
					break;
				case S_FIELD:
					std::cout << S[i*num_cells_y + j] << " ";
					break;
				case P_FIELD:
					std::cout << P[i*num_cells_y + j] << " ";
					break;
				case M_FIELD:
					std::cout << M[i*num_cells_y + j] << " ";
					break;
			}
			
		}
		std::cout << std::endl;
	}
}

void FluidDomain::print_msg(std::string msg)
{
	std::cout << msg << std::endl;
}