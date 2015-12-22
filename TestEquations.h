#ifndef TestEquations_h
#define TestEquations_h
#include <iostream>
#include <string>

#include <cppunit/TestCase.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCaller.h>
#include <cppunit/ui/text/TestRunner.h>
#include "system.h"
#include "sylvester.h"
#include "solver.h"

class EquationsTestCase : public CppUnit::TestCase {
public:
	
	EquationsTestCase() {}
	EquationsTestCase(std::string name) : CppUnit::TestCase(name) {}
	
	void test_calculate_k();
	
	void test_solving();

	void test_change_hidden();
	
	static CppUnit::Test *suite();

};
#endif