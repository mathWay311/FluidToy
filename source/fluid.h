#include <vector>

class FluidDomain
{
public:
	enum field_type {U_FIELD, V_FIELD, S_FIELD};
	enum debug_mode {NONE, MIN, FULL};
	const double overRelaxation = 1.9;
	unsigned int num_cells_x; 
	unsigned int num_cells_y;
	unsigned int num_cells;
	unsigned short density;
	std::vector<double> U;
	std::vector<double> V;
	std::vector<double> newU;
	std::vector<double> newV;
	std::vector<double> P;
	std::vector<double> S;
	std::vector<double> M;
	std::vector<double> newM;

	debug_mode mode;
	unsigned int iteration_counter;

	FluidDomain(unsigned int numx, unsigned int numy, debug_mode mode);

	void integrate(double dt, double gravity);
	void solveIncompressibility(unsigned int numIters, double dt);
	void extrapolate();
	double sampleField(double x, double y, field_type field);
	double avgU(unsigned int i, unsigned int j);
	double avgV(unsigned int i, unsigned int j);
	void advectSmoke(double dt);
	void advectVel(double dt);

	void debug_print_field(field_type type);
	void simulate(double dt, double gravity, unsigned int numIters);

private:
	void print_msg(std::string msg);
};