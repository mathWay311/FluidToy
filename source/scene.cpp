
#include <iostream>
#include "scene.h"

#include "utilities.h"

using enum field_type;

Scene::Scene(int _numX, int _numY)
{
	numX = _numX + 2;
	numY = _numY + 2;
 	fluidDomain = FluidDomain(_numX, _numY, FULL);
 	for (int i = 0; i < numX; i++) 
 	{
 		for (int j = 0; j < numY; j++) 
 		{
				fluidDomain.S[i*numY + j] = 1.0;
				if (i == 0 || i == numX-1 || j == 0) fluidDomain.S[i*numY + j] = 0.0;
				
		}
	}
}

void Scene::add_boundary_inlet(double velocity)
{
	for (int i = 0; i < numX; i++) 
 	{
 		for (int j = 0; j < numY; j++) 
 		{
 			if (i == 1) {
				fluidDomain.U[i*numY + j] = velocity;
				fluidDomain.M[i*numY + j] = 0.0;
			}

 		}
 	}
}

void Scene::simulate()
{
	fluidDomain.simulate(dt, gravity, numIters);
}

void Scene::debug_print()
{
	std::cout << "===\t Поле U\t ===\n" ;
	fluidDomain.debug_print_field(U_FIELD);
	std::cout << "===\t Поле V\t ===\n" ;
	fluidDomain.debug_print_field(V_FIELD);
	std::cout << "===\t Поле S\t ===\n" ;
	fluidDomain.debug_print_field(S_FIELD);
	std::cout << "===\t Поле P\t ===\n" ;
	fluidDomain.debug_print_field(P_FIELD);
	std::cout << "===\t Поле M\t ===\n" ;
	fluidDomain.debug_print_field(M_FIELD);
}