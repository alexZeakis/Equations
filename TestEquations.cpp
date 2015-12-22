#include "TestEquations.h"
#define LIMIT 0.00001

void EquationsTestCase::test_calculate_k() {  

	char * arguments[9];
	strcmp(arguments[0],"./equations");
	strcmp(arguments[1], "-read");
	strcmp(arguments[2], "-i");
	strcmp(arguments[3], "./pols/pol.txt");
	strcmp(arguments[4], "-d1");
	strcmp(arguments[5],"3");
	strcmp(arguments[6],"-d2");
	strcmp(arguments[7], "3");
	strcmp(arguments[8],"-solve");
	sys s1(arguments);
	sylvester syl1(s1);
	double k = syl1.calculate_k();
	assertEquals(-1.0, k);

	strcmp(arguments[3], "./pols/pol4.txt");
	strcmp(arguments[5], "2");
	sys s2(arguments);
	sylvester syl2(s2);
	double k = syl2.calculate_k();
	assert((k > LIMIT && k != -1.0));

}


void EquationsTestCase::test_solving() {
	
	char * arguments[9];
	strcmp(arguments[0], "./equations");
	strcmp(arguments[1], "-read");
	strcmp(arguments[2], "-i");
	strcmp(arguments[3], "./pols/pol5.txt");
	strcmp(arguments[4], "-d1");
	strcmp(arguments[5], "2");
	strcmp(arguments[6], "-d2");
	strcmp(arguments[7], "3");
	strcmp(arguments[8], "-solve");
	sys s1(arguments);
	sylvester syl1(s1);
	solver solve(&syl1, 9, arguments);
	
}

void EquationsTestCase::test_change_hidden(){


}


Test *EquationsTestCase::suite() {
	
	TestSuite *testSuite = new TestSuite("EquationsTestCase");

	
	testSuite->addTest(new TestCaller("test_calculate_k", &EquationsTestCase::test_calculate_k));
	//testSuite->addTest(new TestCaller("test_solving", &EquationsTestCase::test_solving));
	//testSuite->addTest(new TestCaller("test_change_hidden", &EquationsTestCase::test_change_hidden));
	
	return testSuite;
}


int main(int argc, char* argv[]) {
	
	if (argc != 2) {
		std::cout << "usage: tester name_of_class_being_test" << std::endl;
		exit(1);
	}

	TestRunner runner;
	runner.addTest(argv[1], EquationsTestCase::suite());
	runner.run(argc, argv);

	return 0;
}