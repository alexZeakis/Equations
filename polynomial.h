#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <string>

using namespace std;


class polynomial {

	private:
		int** cons;
		int d;

	public:
		polynomial(int d,string pol);
		~polynomial();

		void print();

		int get_d(char var);
		void get_cons(int& dest, int j, int k, char var, int depth);

};

#endif
