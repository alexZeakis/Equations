#ifndef SYLVESTER_H
#define SYLVESTER_H

#include "system.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

class sylvester {

	private:
		int*** smatrix;
		int*** spol;
		int d0,d1, depth;
		double k;

		void print_2d(int **m);

	public:
		sylvester(ssystem& s);
		~sylvester();

		void print_matrix();
		void print_pol(int k);
		int spol_row_sum(int k, int line, int* v);
		double calculate_k();
		void get_spol(MatrixXd& dest, int matrix);

		int getD0(){return d0;};
		int getD1(){return d1;};
		int getDepth(){return depth;};
};
#endif
