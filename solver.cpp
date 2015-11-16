#include <iostream>
#include "solver.h"

solver::solver(sylvester& s, int argc, char *argv[]) {

	k = s.calculate_k();                 /*Claculate k*/
	int b = 7;
	if (strcmp("-solve", argv[argc - 1])){
		b = atoi(argv[argc - 1]);
	}
	/*Case 1 : k < 10^B */
	if (k < pow(10, b) && k > -1){
		cout << endl << "K ~~ " << k << " < Bound: non-singular Μd, standard eigenproblem " << endl;
		c = new companion(s);
		c->solve();
	}
	/*Case 2 : k > 10^B */
	else {
		cout << endl << "K ~~ " << k << " > Bound: ill-conditioned Μd, generalized eigenproblem " << endl;
		l = new  lmatrix(s);
		l->solve();
	}

}

solver::~solver() {
	if(c!=NULL) delete c;
	if(l!=NULL) delete l;
}
