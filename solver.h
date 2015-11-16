#ifndef SOLVER_H
#define SOLVER_H

#include "sylvester.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

class solver {

private:
	

public:
	solver(sylvester& s);
	~solver();

	
};
#endif
