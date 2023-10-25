#include <iostream>
#include "scene.h"
#include <cmath>

Scene::Scene() {}

void Scene::init()	
{
	fluid = FluidDomain(numX, numY, FluidDomain.FULL);
}