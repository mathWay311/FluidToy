#include <iostream>
#include <cmath>

class Scene
{
	public:
	float overRelaxation = 1.9;
	double dt = 1.0 / 60.0;
	unsigned int numIters = 40;

	unsigned int resolution = 100;

	float domainHeight = 1.0;

	int numX;
	int numY;

	float density = 1000.0;
	FluidDomain fluid;

	void init();
	Scene();
}