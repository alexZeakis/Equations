#ifndef SYLVESTER_H
#define SYLVESTER_H

#include "system.h"

using namespace std;

class sylvester {

	private:
		int*** smatrix;
		int*** spol;
		int d0,d1, depth;

		void print_2d(int **m);

	public:
		sylvester(system& s);
		~sylvester();

		void print_matrix();
		void print_pol(int k);
		int spol_row_sum(int k, int line, int* v);
		int getD0(){return d0;};
		int getD1(){return d1;};
		int getDepth(){return depth;};
};
#endif