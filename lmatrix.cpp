#include "lmatrix.h"

lmatrix::lmatrix(sylvester& syl) {
	this->d = syl.getDepth()-1;
	this->m = syl.getD0() + syl.getD1();
	cout << "D is " << d << " and m is " << m << endl;

	l0 = new MatrixXd(d*m, d*m);
	l1 = new MatrixXd(d*m, d*m);
	MatrixXd data(m,m);

/*
	if(d==1) {
//		l1->block(0,0,d*m,d*m) = MatrixXd::Identity(d*m, d*m);
		syl.get_spol(data,d);
		l1->bottomRightCorner(m,m) = data;
		cout << endl << *l1 << endl;
	}
	else {
		l1->block(0,0,d*m,d*m) = MatrixXd::Identity(d*m, d*m);
		syl.get_spol(data,d);
		l1->bottomRightCorner(m,m) = data;
		cout << endl << *l1 << endl;

		l0->topLeftCorner((d-1)*m,m) = MatrixXd::Zero((d-1)*m, m);
		l0->topRightCorner((d-1)*m,(d-1)*m) = (-1)*MatrixXd::Identity((d-1)*m,(d-1)*m);

		for(int i=0; i<d; i++) {
			syl.get_spol(data,i);
			l0->block((d-1)*m, m*i, m, m) = data;
		}
	}
*/
	l1->block(0,0,d*m,d*m) = MatrixXd::Identity(d*m, d*m);
	syl.get_spol(data,d);
	int sign = (d==1)?-1:1;
	l1->bottomRightCorner(m,m) = sign*data;
	cout << endl << *l1 << endl;

	l0->topLeftCorner((d-1)*m,m) = MatrixXd::Zero((d-1)*m, m);
	l0->topRightCorner((d-1)*m,(d-1)*m) = (-1)*MatrixXd::Identity((d-1)*m,(d-1)*m);
	for(int i=0; i<d; i++) {
		syl.get_spol(data,i);
		l0->block((d-1)*m, m*i, m, m) = data;
	}
	cout << endl << *l0 << endl;

}

lmatrix::~lmatrix() {
	delete l0;
	delete l1;
}

void lmatrix::solve() {
//	if(d==1) {
//	}
//	else {

	GeneralizedEigenSolver<MatrixXd> ges;
//	ges.compute(*l0,*l1,true);
	ges.compute(*l1,*l0,true);
/*
	cout << endl << "Eigenvalues are " << endl << ges.eigenvalues() << endl;
	cout << endl << "Eigenvectors are " << endl << ges.eigenvectors() << endl;

	for(int i=0; i<m; i++) {
		if(abs(ges.eigenvalues()[i].real()) >= 0.00001) {
			cout << "y = " << ges.eigenvalues()[i].real() << ", x = " << ges.eigenvectors()(m-2,i).real()/ges.eigenvectors()(m-1,i).real() << endl;
		}
	}
*/
}
