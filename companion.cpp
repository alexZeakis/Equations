#include "companion.h"

companion::companion(sylvester& syl) {
	this->d = syl.getDepth()-1;
	this->m = syl.getD0() + syl.getD1();


//	c = new MatrixXd(d*m, d*m);
	this->c = new MatrixXd(d*m, d*m);
	MatrixXd data(m,m);


	c->topLeftCorner((d-1)*m,m) = MatrixXd::Zero((d-1)*m, m);
	c->topRightCorner((d-1)*m,(d-1)*m) = MatrixXd::Identity((d-1)*m,(d-1)*m);

	syl.get_spol(data, d);
	MatrixXd inv = data.inverse();

	for(int i=0; i<d; i++) {
		syl.get_spol(data,i);
		c->block((d-1)*m, m*i, m, m) = -inv*data;
	}

}

companion::~companion() {
	delete c;
	for(int i=0; i<m; i++)
		solutions[i];
	delete solutions;
}

void companion::solve() {
	EigenSolver<MatrixXd> es(*c);

	solutions = new double*[m];
	for(int i=0; i<m; i++)
		solutions[i] = new double[3];

	for(int i=0; i<m; i++)
		for(int j=0; j<3; j++)
			solutions[i][j] = 0.0;	//initializing solutions

	int bottom=0;	//index of last solution
	for(int i=0; i<m; i++) {
		if(abs(es.eigenvalues()[i].real()) >= 0.00001) {	//solution not close to 0
			int pos = this->find(abs(es.eigenvalues()[i].real()));	//find if this solution has already been found, eg multiplicity more than 1
			if(pos<0) {	//if not
				solutions[bottom][0] = es.eigenvalues()[i].real();	//insert solution-y to solutions
				solutions[bottom][1] = es.eigenvectors()(m-2,i).real()/es.eigenvectors()(m-1,i).real();	//insert solution-x
				pos = bottom;	//convenient storing
				bottom++;	//update index of last
			}
			solutions[pos][2]++;	//increase multiplicity
		}
	}
	this->print();
}

/*function to check if solution has already been found*/
int companion::find(double y) {
	for(int i=0; i<m; i++)
		if( solutions[i][0] == y)
			return i;
	return -1;
}

/*printing function*/
void companion::print() {
	cout << endl << "Roots" << endl << "_______" << endl << endl;
	for(int i=0; i<m; i++) {
		if(solutions[i][0] == 0.0)
			return;
		cout << "y = " << solutions[i][0];
		if(solutions[i][2] > 1.0)
			cout << ", multiplicity = " << solutions[i][2] << endl;
		else
			cout << ", x = " << solutions[i][1] << endl;
	}
}
