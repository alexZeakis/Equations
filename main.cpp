#include <iostream>
#include "system.h"
#include "sylvester.h"
#include "vsylvester.h"
using namespace std;

int main(int argc, char* argv[]) {

	system s(argv);
	s.print();


	sylvester syl(s);
	syl.print_matrix();
	cout << endl;
	syl.print_pol(-1);

	int v[4] = {1,2,3,4};
	vsylvester vsyl(syl,v);
	vsyl.print_matrix();
}
