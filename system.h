#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <string>
#include "polynomial.h"
using namespace std;

class ssystem {

	private:
		int*** sys;
		polynomial** p;
		int col, depth;
		char hidden;

	public:
		ssystem(char* argv[]);
		~ssystem();
		void print();
		int get_d(int pol);
		void get_sys(int** dest, int pol, int skip);
		int getDepth();
};

#endif
