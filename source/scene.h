
#pragma once

#include <iostream>
#include <cmath>
#include "fluid.h"

class Scene
{
	public:
		FluidDomain fluidDomain;
		int numX;
		int numY;
		double dt = 1.0 / 60.0;
		double gravity = 0;
		int numIters = 5;

		Scene(int _numX, int _numY);
		
	public:
		void add_boundary_inlet(double velocity);	
		void simulate();
		void debug_print();
};
