#include "fluid.h"

FluidDomain::FluidDomain(unsigned int numx, unsigned numy)
{
	num_cells_x = numx + 2;
	num_cells_y = numy + 2;
	num_cells = num_cells_y * num_cells_x;
	U.resize(num_cells);
	V.resize(num_cells);
	P.resize(num_cells);
	S.resize(num_cells);
	M.resize(num_cells);
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