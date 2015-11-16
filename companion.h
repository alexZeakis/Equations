#ifndef COMPANION_H
#define COMPANION_H

#include "sylvester.h"

using namespace std;

class companion{

	private:
//		MatrixXd* com;
		int m, d;

	public:
		companion(sylvester& syl);
		~companion();
//		void solve();
};

#endif

