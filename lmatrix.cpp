#include "lmatrix.h"
#include <lapacke.h>

lmatrix::lmatrix(sylvester& syl) {
	this->d = syl.getDepth()-1;
	this->m = syl.getD0() + syl.getD1();
	this->hidden = syl.getHidden();

	l0 = new MatrixXd(d*m, d*m);
	l1 = new MatrixXd(d*m, d*m);
	MatrixXd data(m,m);

	l1->block(0,0,d*m,d*m) = MatrixXd::Identity(d*m, d*m);
	syl.get_spol(data,d);
	int sign = (d==1)?-1:1;
	l1->bottomRightCorner(m,m) = sign*data;

	l0->topLeftCorner((d-1)*m,m) = MatrixXd::Zero((d-1)*m, m);
	l0->topRightCorner((d-1)*m,(d-1)*m) = (-1)*MatrixXd::Identity((d-1)*m,(d-1)*m);
	for(int i=0; i<d; i++) {
		syl.get_spol(data,i);
		l0->block((d-1)*m, m*i, m, m) = data;
	}

}

lmatrix::~lmatrix() {
	delete l0;
	delete l1;

	for(int i=0; i<d*m; i++)
		delete [] solutions[i];
	delete [] solutions;
}

void lmatrix::solve() {
	MatrixXd v(d*m, d*m);
	MatrixXd lambda(d*m, 3);

	int LDA = l0->outerStride();
	int LDB = l1->outerStride();
	int LDV = v.outerStride();
	int INFO= 0;

	double* alphar = lambda.col(0).data();
	double* alphai = lambda.col(1).data();
	double* beta = lambda.col(2).data();

	INFO = LAPACKE_dggev(LAPACK_COL_MAJOR, 'N', 'V', m, l0->data(), LDA, l1->data(), LDB, alphar, alphai, beta, 0, LDV, v.data(), LDV);

	solutions = new double*[m];
	for(int i=0; i<d*m; i++)
		solutions[i] = new double[3];

	for(int i=0; i<d*m; i++)
		for(int j=0; j<3; j++)
			solutions[i][j] = 0.0;	//initializing solutions


//	cout << lambda << endl << endl;
//	cout << v << endl;

	int bottom=0;	//index of last solution
	for(int i=0; i<d*m; i++) {
		if(abs(lambda(i,2)) >= LIMIT) {
			double temp = lambda(i,0)/lambda(i,2);
			if(abs(temp) >= LIMIT) { //hidden.real not close to 0
				int pos = this->find(abs(temp));	//find if this solution has already been found, eg multiplicity more than 1
				if(pos<0) {	//if not
					solutions[bottom][0] = temp; //insert solution-hidden to solutions
					pos = bottom; //convenient storing
					bottom++;	//update index of last
				}
				solutions[pos][1] = i;
				solutions[pos][2]++;	//increase multiplicity
			}
		}
	}

	for(int j=0; j<d*m; j++) {
		if(solutions[j][2] == 1) {
			int i= solutions[j][1];
			if(abs(v(m-1,i)) >= LIMIT && abs(v(m-2,i)) <= LIMIT)
				solutions[j][1] = v(m-2,i)/v(m-1,i);	//insert solution-not-hidden
			else
				solutions[j][2] = 0;
		}
	}
	this->print_solutions();

}

/*function to check if solution has already been found*/
int lmatrix::find(double y) {
	for(int i=0; i<d*m; i++) {
		if( solutions[i][0] == y)
			return i;
	}
	return -1;
}

/*printing solutions function*/
void lmatrix::print_solutions() {
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

void lmatrix::print() {
	cout << "d is " << d << " and m is " << m << endl;

	cout << "L0:" << endl;
	for(int i=0; i<d*m; i++) {
		if(i%m==0) {
			for(int k=0; k<2*d*m; k++)
				cout << "_";
			cout << endl;
		}

		for(int j=0; j<d*m; j+=m) {
			cout << "|" << l0->block(i,j,1,m);
		}
		cout << "|" << endl;
	}

	cout << "L1:" << endl;
	for(int i=0; i<d*m; i++) {
		if(i%m==0) {
			for(int k=0; k<2*d*m; k++)
				cout << "_";
			cout << endl;
		}

		for(int j=0; j<d*m; j+=m) {
			cout << "|" << l1->block(i,j,1,m);
		}
		cout << "|" << endl;
	}

}
