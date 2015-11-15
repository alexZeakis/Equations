#include "companion.h"

companion::companion(sylvester& syl) {
	this->d = syl.getDepth()-1;
	this->m = syl.getD0() + syl.getD1();


//	c = new MatrixXd(d*m, d*m);
	MatrixXd c(d*m, d*m);
	MatrixXd data(m,m);

	c.topLeftCorner((d-1)*m,m) = MatrixXd::Zero((d-1)*m, m);
	c.topRightCorner((d-1)*m,(d-1)*m) = MatrixXd::Identity((d-1)*m,(d-1)*m);

	syl.get_spol(data, d);
//	cout << data << endl;
	MatrixXd inv = data.inverse();
//	cout << inv << endl;
		
	for(int i=0; i<d; i++) {
		syl.get_spol(data,i);
		c.block((d-1)*m, m*i, m, m) = -inv*data;
	}

	cout << endl << c << endl;

	MatrixXd nc = c.transpose();
	cout << endl << nc << endl;

//	EigenSolver<MatrixXd> es(c,true);
	EigenSolver<MatrixXd> es(c);
//	EigenSolver<MatrixXd> es(nc);
//	EigenSolver<MatrixXd> es(nc,true);

	cout << endl << "Solutions-eigenvalues are " << endl << es.eigenvalues() << endl;
	cout << endl << "Eigenvectors are " << endl << es.eigenvectors() << endl;

	for(int i=0; i<m; i++) {
		if(es.eigenvalues()[i].real()!=0) {
//			cout << "y = " << es.eigenvalues()[i].real() << ", x = " << es.eigenvectors()[i][1].real() << endl;
			cout << "y = " << es.eigenvalues()[i].real() << ", x = " << es.eigenvectors()(1,i).real() << endl;
//			cout << es.eigenvectors()(1,i) << endl;
		}
	}
}

companion::~companion() {

}
/*
void companion::solve() {
}
*/
