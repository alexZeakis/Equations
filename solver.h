#ifndef SOLVER_H
#define SOLVER_H

#include "sylvester.h"
#include "companion.h"
#include "lmatrix.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

class solver {

private:
	double k;
	companion* c=NULL;
	lmatrix* l=NULL;


public:
	solver(sylvester& s, int argc, char* argv[]);
	~solver();


};
#endif
