#ifndef COMPANION_H
#define COMPANION_H

#include "sylvester.h"

using namespace std;

class companion{

	private:
		MatrixXd* c;
		double** solutions;
		int m, d;
		char hidden;

		int find(double y);
	public:
		companion(sylvester* syl);
		~companion();
		void solve(int t[] = NULL, char original_hidden = 'y');
		void print_solutions();
		void check_pol_values(sys *s);
		void print();
};

#endif

