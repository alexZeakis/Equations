#include <iostream>
#include "solver.h"
#include "companion.h"


solver::solver(sylvester& s, int argc, char *argv[]) {

	k = s.calculate_k();                 /*Claculate k*/
	int b = 7;
	if (strcmp("-solve", argv[argc - 1])){
		b = atoi(argv[argc - 1]);
	}
	/*Case 1 : k < 10^B */
	if (k < pow(10, b) && k > -1){
		cout << "K ~~ " << k << " < Bound: non-singular	Μd, standard eigenproblem " << endl;
		companion companion_matrix(s);
	}
	/*Case 2 : k > 10^B */
	else {
		cout << "K ~~ " << k << " > Bound: ill-conditioned Μd, generalized eigenproblem " << endl;
		if (s.getDepth() == 1){      /*degree of hidden = 1*/

			int d = s.getDepth() - 1;
			int m = s.getD0() + s.getD1();

			MatrixXd m1(m, m);
			MatrixXd m0(m, m);
			s.get_spol(m0, 0);
			s.get_spol(m1, 1);

			GeneralizedEigenSolver<MatrixXd> ges;
			
			ges.compute(m1, m0);
			cout << "The (complex) numerators of the generalzied eigenvalues are: " << ges.alphas().transpose() << endl;
			cout << "The (real) denominatore of the generalzied eigenvalues are: " << ges.betas().transpose() << endl;
			cout << "The (complex) generalzied eigenvalues are (alphas./beta): " << ges.eigenvalues().transpose() << endl;

			/*for (int i = 0; i<m; i++) {
				if (abs(ges.eigenvalues()[i].real()) >= 0.00001) {
					cout << "y = " << ges.eigenvalues()[i].real() << ", x = " << ges.eigenvectors()(1, i).real() / ges.eigenvectors()(2, i).real() << endl;
				}
			}*/
		}
		else{
			cout << "Generalized eigenproblem with d >= 2 " << endl;
			//lmatrix l_matrix(s);
		}
	}

}

solver::~solver() {

	cout << "Solver destructor" << endl;

}