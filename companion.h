#ifndef COMPANION_H
#define COMPANION_H

#include "sylvester.h"

using namespace std;

class companion{

	private:
		MatrixXd* c;
		double** solutions;
		int m, d;

		void print();
		int find(double y);
	public:
		companion(sylvester& syl);
		~companion();
		void solve();
};

#endif

