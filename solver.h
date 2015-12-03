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
	int b;
	companion* c = NULL;
	lmatrix* l = NULL;
	sylvester* syl;
	int solve(int t[] = NULL, char original_hidden = 'y');


public:
	solver(sylvester* s, int argc, char* argv[]);
	~solver();

	int change_hidden();
	


};
#endif