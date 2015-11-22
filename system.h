#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <string>
#include "polynomial.h"
using namespace std;

class sys {

	private:
		int*** system_matrix;
		polynomial** p;
		int col, depth;
		char hidden;

	public:
		sys(char* argv[]);
		~sys();
		void print();
		int get_d(int pol);
		void get_sys(int** dest, int pol, int skip);
		int getDepth();
		char getHidden() {return hidden;};
};

#endif
