#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <string>

using namespace std;


class polynomial {

	private:
		double** cons;
		int d;

	public:
		polynomial(int d,string pol);
		~polynomial();

		void print();
		void printGenerate();

		int get_d(char var);
		void get_cons(double& dest, int j, int k, char var, int depth);
		double calculate_value(double x, double y);

};

#endif
