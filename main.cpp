#include <iostream>
#include "system.h"
#include "sylvester.h"
#include "vsylvester.h"
#include "companion.h"
#include "solver.h"
#include <Eigen/Dense>

using namespace std;

int main(int argc, char* argv[]) {

	cout << endl << "Equations" << endl << "__________________" << endl;
	sys s(argv);                            /*Make system */

	//s.print();

	sylvester syl(s);                       /*Make sylvester matrix */

	//	syl.print_matrix();
	//	cout << endl;
	//	syl.print_pol(-1);

	syl.print_pol(syl.getDepth() - 1);     /*Print Md */

	//	cout << "K ~~ " << syl.calculate_k() << endl;

	solver solve(&syl, argc, argv);        /* Solve */
//	solve.change_hidden();                 /* Attempt to change hidden variable, possibly solve again */

}
