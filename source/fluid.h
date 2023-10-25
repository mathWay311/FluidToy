#include <vector>
class FluidDomain
{
public:
	unsigned int num_cells_x; 
	unsigned int num_cells_y;
	unsigned int num_cells;
	unsigned short density;
	std::vector<double> U;
	std::vector<double> V;
	std::vector<double> P;
	std::vector<double> S;
	std::vector<double> M;

	FluidDomain(unsigned int numx, unsigned numy);

	void integrate(double dt, double gravity);


};