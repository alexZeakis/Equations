#include <iostream>
#include "sylvester.h"

sylvester::sylvester(sys& s) {
	this->d0 = s.get_d(0);
	this->d1 = s.get_d(1);
	this->depth = s.getDepth();
	this->hidden = s.getHidden();

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

	int test[3]= {-1,-1,-1};
	MatrixXd testm(d0+d1,d0+d1);
	MatrixXd My(d0+d1, d0+d1);
	My = MatrixXd::Zero(d0+d1, d0+d1);
	for(int i=0; i<3; i++) {
		int temp;
		do {
			temp = rand() % 100;
		} while(temp==test[0] || temp==test[1] || temp==test[2]);	//different tests
		test[i] = temp;
		for(int j=0; j<depth; j++) {
			this->get_spol(testm,j);	//useless as an accessor, but good as a converter from int** to MatrixXd
			My = My + testm*pow(temp,j);
		}
		if(My.determinant()!= 0)	//no need for further tests
			break;
		else if(i==2) {	//all determinants were 0 and this is the 3rd test
			cout << "Determinants were all 0" << endl;
			exit(2);
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

void sylvester::get_spol(MatrixXd& dest, int matrix) {
	for(int i=0; i<d0+d1; i++)
		for(int j=0; j<d0+d1; j++)
			dest(i,j) = spol[matrix][i][j];
}

double sylvester::calculate_k() {
	MatrixXd data(d0+d1, d0+d1);
	this->get_spol(data, depth - 1);          /*Write Md from the sylvester polyonym matrix into matrix data*/

	JacobiSVD<MatrixXd> svd(data, ComputeFullU);

	double min = svd.singularValues()[0];
	double max = svd.singularValues()[0];

	for (int i = 0; i< d0 + d1; i++) {          /*Find σmax and σmin*/
		double temp = svd.singularValues()[i];
		if(temp < min)
			min = temp;
		if(temp > max)
			max = temp;

	}

	return k = (min < LIMIT)?-1:max/min; 	// min is close to 0 => k tends to infinity
}
