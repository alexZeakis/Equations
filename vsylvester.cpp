#include "vsylvester.h"
#include <stdlib.h>

vsylvester::vsylvester(sylvester& syl, int* v) {
	this->d0 = syl.getD0();
	this->d1 = syl.getD1();
	if(d0+d1!= sizeof(v)) {
		cout << "Wrong Dimension of vector v. It should be " << d0+d1 << " but it is " << sizeof(v) << endl;
		exit(2);
	}

	this->depth = syl.getDepth();

	this->smatrix = new int*[d0+d1];
	for(int i=0; i<d0+d1; i++) {
		this->smatrix[i] = new int[depth];
	}

	for(int i=0; i<d0+d1; i++) {
		for(int j=0; j<depth; j++) {
			smatrix[i][j] = syl.spol_row_sum(j,i,v);
		}
	}

}

vsylvester::~vsylvester() {
	for(int i=0; i<d0+d1; i++)
		delete smatrix[i];
	delete smatrix;
}

void vsylvester::print_matrix() {
	for(int i=0; i<d0+d1; i++) {
		cout << "[" ;
		for(int j=0; j<depth-1; j++) {
			cout << smatrix[i][j] << ",";
		}
		cout << smatrix[i][depth-1] <<"]" << endl;
	}
}
