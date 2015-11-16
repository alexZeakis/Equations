#ifndef VSYLVESTER_H
#define VSYLVESTER_H

#include "sylvester.h"

using namespace std;

class vsylvester {

	private:
		int** smatrix;
		int d0,d1, depth;

	public:
		vsylvester(sylvester& syl,int* v);
		~vsylvester();

		void print_matrix();

};
#endif
