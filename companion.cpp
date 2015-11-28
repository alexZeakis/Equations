#include "companion.h"

companion::companion(sylvester& syl) {
	this->d = syl.getDepth()-1;
	this->m = syl.getD0() + syl.getD1();
	this->hidden = syl.getHidden();

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
	for(int i=0; i<d*m; i++)
		delete [] solutions[i];
	delete [] solutions;
}

void companion::solve() {
	EigenSolver<MatrixXd> es(*c);

	solutions = new double*[d*m];
	for(int i=0; i<d*m; i++)
		solutions[i] = new double[3];

	for(int i=0; i<d*m; i++)
		for(int j=0; j<3; j++)
			solutions[i][j] = 0.0;	//initializing solutions

	int bottom=0;	//index of last solution
	for(int i=0; i<d*m; i++) {
		complex<double> temp = es.eigenvalues()[i];
		if(abs(temp.real()) >= LIMIT && abs(temp.imag()) <= LIMIT) {	//hidden.real not close to 0, hidden.imag realy small
			int pos = this->find(temp.real());	//find if this solution has already been found, eg multiplicity more than 1
			if(pos<0) {	//if not
				solutions[bottom][0] = temp.real();	//insert solution-hidden to solutions
				pos = bottom;	//convenient storing
				bottom++;	//update index of last
			}
			solutions[pos][1] = i; // convenient index stroring
			solutions[pos][2]++;	//increase multiplicity
		}
	}

	for(int j=0; j<d*m; j++) {
		if(solutions[j][2] == 1) {
			int i= solutions[j][1];
			if(abs(es.eigenvectors()(m-1,i).real()) >= LIMIT && abs(es.eigenvectors()(m-1,i).imag()) <= LIMIT)  //not-hidden.real not close to 0, not-hidden.imag realy small
				solutions[j][1] = es.eigenvectors()(m-2,i).real()/es.eigenvectors()(m-1,i).real();	//insert solution-not-hidden
			else
				solutions[j][2] = 0;
		}
	}
	this->print_solutions();
}

/*function to check if solution has already been found*/
int companion::find(double y) {
	for(int i=0; i<d*m; i++)
		if( solutions[i][0] == y)
			return i;
	return -1;
}

/*printing solutions function*/
void companion::print_solutions() {
	cout << endl << "Roots" << endl << "_______" << endl << endl;
	for(int i=0; i<d*m; i++) {
		if(solutions[i][2] > 0) {
			cout << hidden << " = " << solutions[i][0];
			if(solutions[i][2] > 1.0)
				cout << ", multiplicity = " << solutions[i][2] << endl;
			else
				cout << ", " << ((hidden=='x')?"y":"x") << " = " << solutions[i][1] << endl;
		}
	}
}

void companion::print() {
	cout << "Companion Matrix: " << endl;
	cout << *c << endl;
}
