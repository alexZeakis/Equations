#include <iostream>
#include "sylvester.h"



sylvester::sylvester(sys& s, sylvester* syl, int t[]) {

	if (syl == NULL){                   /* CASE 1 : Making sylvester matrix from system object */
		this->d0 = s.get_d(0);
		this->d1 = s.get_d(1);
		this->depth = s.getDepth();
		this->hidden = s.getHidden();
		this->origin_system = &s;

		smatrix = new double**[d0 + d1];
		for (int i = 0; i < d0 + d1; i++) {
			smatrix[i] = new double*[d0 + d1];
			for (int j = 0; j < d0 + d1; j++) {
				smatrix[i][j] = new double[depth];
			}
		}

		for (int i = 0; i < d0 + d1; i++) {
			for (int j = 0; j < d0 + d1; j++) {
				for (int k = 0; k < depth; k++) {
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


		/* Fill spol, representing sylvester polyonym from smatrix */

		spol = new double**[depth];
		for (int i = 0; i < depth; i++) {
			spol[i] = new double*[d0 + d1];
			for (int j = 0; j < d0 + d1; j++) {
				spol[i][j] = new double[d0 + d1];
			}
		}

		for (int i = 0; i < depth; i++) {
			for (int j = 0; j < d0 + d1; j++) {
				for (int k = 0; k < d0 + d1; k++) {
					spol[i][j][k] = smatrix[j][k][i];
				}
			}
		}

		/*_______ Bonus 5 _____________________________*/

		int test[3] = { -1, -1, -1 };
		MatrixXd testm(d0 + d1, d0 + d1);
		MatrixXd My(d0 + d1, d0 + d1);
		My = MatrixXd::Zero(d0 + d1, d0 + d1);
		for (int i = 0; i<3; i++) {
			int temp;
			do {
				temp = rand() % 100;
			} while (temp == test[0] || temp == test[1] || temp == test[2]);	/* 3 different tests */
			test[i] = temp;
			for (int j = 0; j<depth; j++) {
				this->get_spol(testm, j);
				My = My + testm*pow(temp, j);
			}
			if (My.determinant() != 0)	/* Found non-zero, stop testing */
				break;
			else if (i == 2) {
				cout << "Determinants were all 0" << endl;
				exit(2);
			}
		}

		/*__________________________________________*/
	}
	else{                                   /* CASE 2 : Making a new sylvester matrix from another one */
		                                    /*          with a change of hidden variable               */
		this->d0 = syl->getD0();
		this->d1 = syl->getD1();
		this->depth = syl->getDepth();
		this->hidden = 'z';
		this->origin_system = syl->getSystem();

		smatrix = new double**[d0 + d1];
		for (int i = 0; i < d0 + d1; i++) {
			smatrix[i] = new double*[d0 + d1];
			for (int j = 0; j < d0 + d1; j++) {
				smatrix[i][j] = new double[depth];
			}
		}

		spol = new double**[depth];
		for (int i = 0; i < depth; i++) {
			spol[i] = new double*[d0 + d1];
			for (int j = 0; j < d0 + d1; j++) {
				spol[i][j] = new double[d0 + d1];
			}
		}

		/* Initialize spol matrix to 0 */
		for (int i = 0; i < depth; i++) {
			for (int j = 0; j < d0 + d1; j++) {
				for (int k = 0; k < d0 + d1; k++) {
					spol[i][j][k] = 0;
				}
			}
		}

		//cout << endl << "t1 = " << t[0] << " t2 = " << t[1] << " t3 = " << t[2] << " t4 = " << t[3] << endl;

		MatrixXd *temp1 = new MatrixXd[depth];
		MatrixXd *temp2 = new MatrixXd[depth];

		/* Initialize temp1[depth][d0 + d1][d0 + d1], temp2[depth][d0 + d1][d0 + d1] to zero*/
		for (int i = 0; i < depth; i++) {
			temp1[i] =  MatrixXd::Zero(d0 + d1, d0 + d1);
			temp2[i] =  MatrixXd::Zero(d0 + d1, d0 + d1);
		}

		for (int i = 0; i < depth; i++) {     /* For every Mi, calculate temp1 (its contribution to every Mi')*/

			MatrixXd& ptr = temp1[0];         /* Copy original Mi into temp1[0]*/
			syl->get_spol(ptr, i);

			for (int j = 0; j < i; j++){      /* Calculate Mi*(t1z + t2)^j and keep it in the appropriate index of temp1 */

				for (int l = 0; l <= j; l++){  /*Use temp2 to save the sum of multiplications */

					MatrixXd& t1 = temp1[l];
					MatrixXd& t2 = temp2[l];
					MatrixXd& t3 = temp2[l + 1];
					t2 += (t1 * t[1]);           /*temp2[l] += temp1[l]*t2 */
					t3 += (t1 * t[0]);           /*temp2[l + 1] += temp1[l]*t1*/

				}

				for (int l = 0; l <= j + 1; l++) { /*Copy the sums saved into temp2 to temp1, return temp2 to zero*/

					MatrixXd& t1 = temp1[l];
					MatrixXd& t2 = temp2[l];
					t1 = t2;
					t2.setZero();
				}
			}
			for (int j = 0; j < depth - 1 - i; j++){

				/* Calculate Mi*(t3z + t4)^j and keep it in the appropriate index of temp1 */

				for (int l = 0; l <= j; l++){

					MatrixXd& t1 = temp1[l];
					MatrixXd& t2 = temp2[l];
					MatrixXd& t3 = temp2[l + 1];
					t2 += (t1 * t[3]);
					t3 += (t1 * t[2]);
				}

				for (int l = 0; l <= j + 1; l++) {

					MatrixXd& t1 = temp1[l];
					MatrixXd& t2 = temp2[l];
					t1 = t2;
					t2.setZero();
				}
			}

			/* Add temp1 (the contribution of Mi to every Mi') to the new spol matrix */
			for (int j = 0; j < depth; j++) {
				for (int k = 0; k < d0 + d1; k++) {
					for (int l = 0; l < d0 + d1; l++) {
						spol[j][k][l] += temp1[j](k, l);
					}
				}
			}

		}

		delete [] temp1;
		delete [] temp2;

		/* Fill smatrix from the sylvester polyonym matrix */
		for (int i = 0; i < depth; i++) {
			for (int j = 0; j < d0 + d1; j++) {
				for (int k = 0; k < d0 + d1; k++) {
					smatrix[j][k][i] = spol[i][j][k];
				}
			}
		}

	}

}

sylvester::~sylvester() {


	for (int i = 0; i < d0 + d1; i++) {
		for (int j = 0; j < d0 + d1; j++) {
			delete[] smatrix[i][j];
		}
		delete[] smatrix[i];
	}
	delete[] smatrix;

	for (int i = 0; i < depth; i++) {
		for (int j = 0; j < d0 + d1; j++) {
			delete[] spol[i][j];
		}
		delete[] spol[i];
	}
	delete[] spol;

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
void sylvester::print_2d(double **matrix) {

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
				cout << hidden;
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

/* Multiply row with vector, return sum of row */
int sylvester::spol_row_sum(int k, int line, int*v) {
	int sum = 0;
	for (int j = 0; j<d0 + d1; j++) {
		sum += spol[k][line][j] * v[j];
	}
	return sum;
}

/* Write Md, d=matrix into matrixXd dest */
void sylvester::get_spol(MatrixXd& dest, int matrix) {
	for (int i = 0; i<d0 + d1; i++)
		for (int j = 0; j<d0 + d1; j++)
			dest(i, j) = spol[matrix][i][j];
}


double sylvester::calculate_k() {
	MatrixXd data(d0 + d1, d0 + d1);
	this->get_spol(data, depth - 1);          /*Write Md from the sylvester polyonym matrix into matrix data*/

	JacobiSVD<MatrixXd> svd(data, ComputeFullU);

	double min = svd.singularValues()[0];
	double max = svd.singularValues()[0];

	for (int i = 0; i< d0 + d1; i++) {          /*Find σmax and σmin*/
		double temp = svd.singularValues()[i];
		if (temp < min)
			min = temp;
		if (temp > max)
			max = temp;

	}

	return k = (min < LIMIT) ? -1 : max / min; 	/*min is close to 0 => k tends to infinity, k = -1*/
}
