#ifndef LMATRIX_H
#define LMATRIX_H

#include "sylvester.h"

using namespace std;

class lmatrix{

	private:
		MatrixXd* l0;
		MatrixXd* l1;
		double ** solutions;
		int m, d;

		void print();
		int find(double y);

	public:
		lmatrix(sylvester& syl);
		~lmatrix();
		void solve();
};

#endif
