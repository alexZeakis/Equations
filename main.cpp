#include <iostream>
#include "system.h"
#include "sylvester.h"
#include "vsylvester.h"
#include <Eigen/Dense>

using namespace std;

int main(int argc, char* argv[]) {

	sys s(argv);
//	s.print();

	sylvester syl(s);
//	syl.print_matrix();
//	cout << endl;
//	syl.print_pol(-1);
	syl.print_pol(syl.getDepth()-1);

	cout << "K is " << syl.calculate_k() << endl;
/*
	int v[4] = {1,2,3,4};
	vsylvester vsyl(syl,v);
	vsyl.print_matrix();
*/
}
