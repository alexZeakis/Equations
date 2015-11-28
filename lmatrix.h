#ifndef LMATRIX_H
#define LMATRIX_H

#include "sylvester.h"

using namespace std;

class lmatrix{

	private:
		MatrixXd* l0;
		MatrixXd* l1;
		double** solutions;
		int m, d;
		char hidden;

		void print_solutions();
		int find(double y);

	public:
		lmatrix(sylvester& syl);
		~lmatrix();
		void solve();
		void print();
};

#endif

