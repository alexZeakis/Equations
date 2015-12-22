#include <iostream>
#include "solver.h"
#define T_RANGE 100;

solver::solver(sylvester* s, int argc, char *argv[]) {

	syl = s;
	k = syl->calculate_k();                 /*Calculate k*/
	b = 7;
	
	if (strcmp("-solve", argv[argc - 1])){  /*If last user argument isn't "-solve", user has given b*/
		b = atoi(argv[argc - 1]);
	}

	this->solve();
}

int solver::solve(int t[], char original_hidden){
		
	/*Case 1 : k < 10^B */
	if (k < pow(10, b) && k > -1){
		cout << endl << "K ~~ " << k << " < Bound: non-singular Μd, standard eigenproblem " << endl;
		c = new companion(syl);
		if (t != NULL)
			c->solve(t, original_hidden);   /* Extra info needed when solving a changed variable problem */
		else
			c->solve();
		c->print_solutions();
		c->check_pol_values(syl->getSystem());
	}
	/*Case 2 : k > 10^B */
	else {
		cout << endl << "K ~~ " << k << " > Bound: ill-conditioned Μd, generalized eigenproblem " << endl;
		l = new  lmatrix(syl);
		if (t != NULL)
			l->solve(t, original_hidden);    /* Extra info needed when solving a changed variable problem */
		else
			l->solve();
		l->print_solutions();
		l->check_pol_values(syl->getSystem());
	}

	return 0;

}

int solver::change_hidden(){

	int attempt = 0;

	while (attempt < 4){        /*Up to 4 attempts to change the hidden variable */

		sys *null_sys = NULL;
		int t[4];
		for (int i = 0; i < 4; i++){           /* Pick 4 random ti */
			t[i] = rand() % T_RANGE + 1;
		}
		sylvester new_syl(*null_sys, syl, t);  /* Create new sylvester, based on the existing one */

//		syl->print_matrix();
//		new_syl.print_matrix();
//		cout << endl;
//		syl->print_pol(-1);
//		new_syl.print_pol(-1);

		double new_k = new_syl.calculate_k();
		cout << endl << "Attempt to change hidden variable resulted in new K ~~  " << new_k << endl;
		
		if (new_k != -1 && (new_k < k || k == -1)){  /* Determine if k is better */

			sylvester *old_syl = syl;        /* Save old sylvester pointer */
			char original_hidden = old_syl->getHidden();
			syl = &new_syl;
			k = new_k;
			if (c != NULL){
				delete c;
				c = NULL;
			}
			if (l != NULL){
				delete l;
				l = NULL;
			}
			this->solve(t, original_hidden);                 /* Solve with new sylvester */
			syl = old_syl;                                   /* Return to original sylvester */
			break;
		}
		attempt++;
	}

	return 0;

}

solver::~solver() {
	if (c != NULL) delete c;
	if (l != NULL) delete l;
}


void solver::print() {
	if (c != NULL) c->print();
	if (l != NULL) l->print();
}