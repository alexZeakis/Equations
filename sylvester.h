#ifndef SYLVESTER_H
#define SYLVESTER_H

#define LIMIT 0.00000000000000000000001
#define TESTLIMIT 0.0000000001
//#define LIMIT 0.000000001

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
		char hidden;
		sys *origin_system;

		void print_2d(int **m);

	public:
		sylvester(sys& s, sylvester* syl = NULL, int t[] = NULL);
		~sylvester();

		void print_matrix();
		void print_pol(int k);
		int spol_row_sum(int k, int line, int* v);
		double calculate_k();
		void get_spol(MatrixXd& dest, int matrix);

		int getD0(){ return d0; };
		int getD1(){ return d1; };
		int getDepth(){ return depth; };
		char getHidden(){ return hidden; };
		sys *getSystem(){ return origin_system; };

};
#endif
