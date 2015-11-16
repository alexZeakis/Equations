#include <iostream>
#include "sylvester.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

sylvester::sylvester(sys& s) {
	this->d0 = s.get_d(0);
	this->d1 = s.get_d(1);
	this->depth = s.getDepth();

	smatrix = new int**[d0+d1];
	for(int i=0; i<d0+d1; i++) {
		smatrix[i] = new int*[d0+d1];
		for(int j=0; j<d0+d1; j++) {
			smatrix[i][j] = new int[depth];
		}
	}

	for(int i=0; i<d0+d1; i++) {
		for(int j=0; j<d0+d1; j++) {
			for(int k=0; k < depth; k++) {
				smatrix[i][j][k] = 0;   /* intializing smatrix of size [d0+d1][d0+d1][max degree of hidden var]*/
			}
		}
	}

	/* Fill smatrix by copying parts of the system matrix in the correct positions */
	for (int i = 0; i < d1; i++){
		s.get_sys(smatrix[i], 0, i);
	}
	for (int i = 0; i < d0; i++){
		s.get_sys(smatrix[d1 + i], 1, i);
	}


	/*__________________________________________*/

	spol = new int**[depth];
	for(int i=0; i < depth; i++) {
		spol[i] = new int*[d0+d1];
		for(int j=0; j<d0+d1; j++) {
			spol[i][j] = new int[d0+d1];
		}
	}

	for (int i = 0; i < depth; i++) {
		for (int j = 0; j < d0 + d1; j++) {
			for (int k = 0; k < d0 + d1; k++) {
				spol[i][j][k] = smatrix[j][k][i];
			}
		}
	}

}

sylvester::~sylvester() {


	for (int i = 0; i < d0 + d1; i++) {
		for (int j = 0; j < d0 + d1; j++) {
				delete smatrix[i][j];			
		}
		delete smatrix[i];
	}
	delete smatrix;

	for (int i = 0; i < depth; i++) {
		for (int j = 0; j < d0 + d1; j++) {
			delete spol[i][j];
		}
		delete spol[i];
	}
	delete spol;

}

/*Print Sylvester matrix */
void sylvester::print_matrix() {

	cout << endl;
	for (int i = 0; i < d0 + d1; i++) {
		for (int j = 0; j < d0 + d1; j++) {
			cout << "[";
			for (int k = 0; k < depth - 1; k++) {
				cout << smatrix[i][j][k] << ",";
			}
			cout << smatrix[i][j][depth - 1] << "]\t";
		}
		cout << endl;
	}
	cout << endl;
}

/*Print a 2D matrix and dont change line after printg the last row */
void sylvester::print_2d(int **matrix) {
	
	cout << endl;
	for (int i = 0; i < d0 + d1; i++) {
		for (int j = 0; j < d0 + d1; j++) {
			cout << matrix[i][j] << "\t";
		}
		if (i != d0 + d1 - 1){
			cout << endl;
		}
	}

}

/* Print the matrix A(k) of the Sylvester polyonym, print the entire Sylvester polyonym for -1 */
void sylvester::print_pol(int k) {

	if (k == -1){

		for (int i = 0; i < depth; i++){
			this->print_2d(spol[i]);
			cout << "\t";
			if (i > 0)
				cout << "y";
			if (i > 1)
				cout << "^" << i;
			if (i < depth - 1){
				cout << "\t+"; 
			}
			cout << endl;
		}
	}
	else{
		if (k > depth){
			cout << "The Sylvester polyonym degree is less than " << k << endl;
			return;
		}
		this->print_2d(spol[k]);
		cout << endl;
	}

	return;


}

/* R*/
int sylvester::spol_row_sum(int k, int line, int*v) {
	int sum=0;
	for(int j=0; j<d0+d1; j++) {
		sum+= spol[k][line][j]*v[j];
	}
	return sum;
}


double sylvester::calculate_k() {
	double *data = new double[(d0+d1)*(d0+d1)];
	for(int i=0; i<d0+d1; i++)
		for(int j=0; j<d0+d1; j++)
			data[(i*(d0+d1))+j] = spol[depth-1][i][j];

	Eigen::Map<Eigen::MatrixXd> m(data, d0+d1, d0+d1);

	Eigen::EigenSolver<Eigen::MatrixXd> es(m);

	cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;

	double min=es.eigenvalues()[0].real();
	double max=es.eigenvalues()[0].real();
	for(int i=0; i< d0+d1; i++) {
		if(es.eigenvalues()[i].real() < min)
			min = es.eigenvalues()[i].real();
		if(es.eigenvalues()[i].real() > max)
			max = es.eigenvalues()[i].real();
	}

	cout << "Min is " << min << endl;
	cout << "Max is " << max << endl;

	delete data;
	return this->k= max/min;
}
