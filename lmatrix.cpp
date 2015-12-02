#include "lmatrix.h"
#include <lapacke.h>

lmatrix::lmatrix(sylvester* syl) {
	this->d = syl->getDepth()-1;
	this->m = syl->getD0() + syl->getD1();
	this->hidden = syl->getHidden();

	l0 = new MatrixXd(d*m, d*m);
	l1 = new MatrixXd(d*m, d*m);
	MatrixXd data(m,m);

	l1->block(0,0,d*m,d*m) = MatrixXd::Identity(d*m, d*m);
	syl->get_spol(data,d);
	int sign = (d==1)?-1:1;
	l1->bottomRightCorner(m,m) = sign*data;

	l0->topLeftCorner((d-1)*m,m) = MatrixXd::Zero((d-1)*m, m);
	l0->topRightCorner((d-1)*m,(d-1)*m) = (-1)*MatrixXd::Identity((d-1)*m,(d-1)*m);
	for(int i=0; i<d; i++) {
		syl->get_spol(data,i);
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

void lmatrix::solve(int t[]) {
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
			solutions[i][j] = 0.0;	/*intialize solutions */


//	cout << lambda << endl << endl;
//	cout << v << endl;

	int bottom=0;	                 /* holds index of last solution */
	for(int i=0; i<d*m; i++) {
		if(abs(lambda(i,2)) >= LIMIT) {
			double temp = lambda(i,0)/lambda(i,2);
			if(abs(temp) >= LIMIT) {  /*hidden.real not close to 0 */

				/* check if this solution has already been found, eg multiplicity more than 1*/
				int pos = this->find(abs(temp));	
				if(pos<0) {	                            /*this is a new solution*/
					solutions[bottom][0] = temp;        /*insert to solutions */
					pos = bottom;                       /*temporarily save bottom */
					bottom++;	                        /*update index of last solution */
				}
				solutions[pos][1] = i;
				solutions[pos][2]++;	                /*increase multiplicity for solution */
			}
		}
	}

	for(int j=0; j<d*m; j++) {
		if(solutions[j][2] == 1) {
			int i= solutions[j][1];
			if(abs(v(m-1,i)) >= LIMIT && abs(v(m-2,i)) <= LIMIT)
				solutions[j][1] = v(m-2,i)/v(m-1,i);	/*insert solution-not-hidden */
			else
				solutions[j][2] = 0;
		}
	}

	if (t != NULL){        /* If hidden variable is z, due to change of variable */
		for (int j = 0; j < d*m; j++){

			solutions[j][0] = (t[0] * solutions[j][0] + t[1]) / (t[2] * solutions[j][0] + t[3]);
		}
		this->hidden = 'y';
	}
	
}

/* Checks if solution has already been found */
int lmatrix::find(double y) {
	for(int i=0; i<d*m; i++) {
		if( solutions[i][0] == y)
			return i;
	}
	return -1;
}

/* Print solutions */
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

/* Print matrices L0, L1 */
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

void lmatrix::check_pol_values(sys *s){

	polynomial *p1 = s->get_pol(0);
	polynomial *p2 = s->get_pol(1);
	for (int i = 0; i < d*m; i++){
		if (solutions[i][2] == 1) {
			double x = (hidden == 'x') ? solutions[i][0] : solutions[i][1];
			double y = (hidden == 'x') ? solutions[i][1] : solutions[i][0];

			if (p1->calculate_value(x, y) < LIMIT && p2->calculate_value(x, y) < LIMIT)
				cout << endl << "Roots x = " << x << ", y = " << y << " is correct when tested on the polynomials. " << endl;
			else
				cout << endl << "Roots x = " << x << ", y = " << y << " is not correct when tested on the polynomials. " << endl;

			//cout << "Check : p1 = " << p1->calculate_value(x, y) << " and p2 = " << p2->calculate_value(x, y) << endl;
		}
	}
	
}
