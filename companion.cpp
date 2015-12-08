#include "companion.h"

companion::companion(sylvester* syl) {
	this->d = syl->getDepth()-1;
	this->m = syl->getD0() + syl->getD1();
	this->hidden = syl->getHidden();
	this->c = new MatrixXd(d*m, d*m);
	MatrixXd data(m,m);

	/* Set up zero and indentity blocks */
	c->topLeftCorner((d-1)*m,m) = MatrixXd::Zero((d-1)*m, m);
	c->topRightCorner((d-1)*m,(d-1)*m) = MatrixXd::Identity((d-1)*m,(d-1)*m);

	syl->get_spol(data, d);
	MatrixXd inv = data.inverse(); /* inv = Md^-1 */

	/* Fill in bottom row with -Md^-1 * Mi */
	for(int i=0; i<d; i++) {
		syl->get_spol(data,i);
		c->block((d-1)*m, m*i, m, m) = -inv*data;
	}

}

companion::~companion() {
	delete c;
	for(int i=0; i<d*m; i++)
		delete [] solutions[i];
	delete [] solutions;
}

void companion::solve(int t[], char original_hidden) {
	EigenSolver<MatrixXd> es(*c);

	solutions = new double*[d*m];
	for(int i=0; i<d*m; i++)
		solutions[i] = new double[3];

	for(int i=0; i<d*m; i++)
		for(int j=0; j<3; j++)
			solutions[i][j] = 0.0;	/* initialize solutions */

	int bottom=0;	                /* holds index of last solution */
	for(int i=0; i<d*m; i++) {
		complex<double> temp = es.eigenvalues()[i];
		if(abs(temp.real()) >= LIMIT && abs(temp.imag()) <= LIMIT) { /*hidden.real not close to 0, hidden.imag really small*/
			
			/* check if this solution has already been found, eg multiplicity more than 1*/
			int pos = this->find(temp.real());	
			if (pos<0) {	                         /*this is a new solution*/
				solutions[bottom][0] = temp.real();	 /*insert to solutions */
				pos = bottom;                        /*temporarily save bottom */
				bottom++;	                         /*update index of last solution */
			}
			solutions[pos][1] = i; 
			solutions[pos][2]++;	                /*increase multiplicity for solution */
		}
	}

	for(int j=0; j<d*m; j++) {
		if(solutions[j][2] == 1) {
			int i= solutions[j][1];
			if(abs(es.eigenvectors()(m-1,i).real()) >= LIMIT && abs(es.eigenvectors()(m-1,i).imag()) <= LIMIT)  
				
				/* If not-hidden.real not close to 0, not-hidden.imag realy small, insert solution not-hidden */
				solutions[j][1] = es.eigenvectors()(m-2,i).real()/es.eigenvectors()(m-1,i).real();	
			else
				solutions[j][2] = 0;
		}
	}

	if (t != NULL){        /* If hidden variable is z, due to change of variable */
		for (int j = 0; j < d*m; j++){
			//if (solutions[j][2] >= 1)
				solutions[j][0] = (t[0] * solutions[j][0] + t[1]) / (t[2] * solutions[j][0] + t[3]);
		}
		this->hidden = original_hidden;
	}
	
}

/* Checks if solution has already been found */
int companion::find(double y) {
	for(int i=0; i<d*m; i++)
		if( solutions[i][0] == y)
			return i;
	return -1;
}

/* Print solutions */
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

/* Print companion matrix */
void companion::print() {
	cout << "Companion Matrix: " << endl;
	cout << *c << endl;
}

void companion::check_pol_values(sys *s){

	polynomial *p1 = s->get_pol(0);
	polynomial *p2 = s->get_pol(1);
	for (int i = 0; i < d*m; i++){
		if (solutions[i][2] == 1) {
			double x = (hidden == 'x') ? solutions[i][0] : solutions[i][1];
			double y = (hidden == 'x') ? solutions[i][1] : solutions[i][0];

			if (p1->calculate_value(x, y) < TESTLIMIT && p2->calculate_value(x, y) < TESTLIMIT)
				cout << endl << "Roots x = " << x << ", y = " << y << " is correct when tested on the polynomials. " << endl;
			else
				cout << endl << "Roots x = " << x << ", y = " << y << " is not correct when tested on the polynomials. " << endl;

			cout << "Check : p1 = " << p1->calculate_value(x, y) << " and p2 = " << p2->calculate_value(x, y) << endl;
		}
	}

}
