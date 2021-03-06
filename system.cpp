#include <stdlib.h>
#include "system.h"
#include <fstream>
#include <string>
#include <string.h>
#include <time.h>
#include <unistd.h>

sys::sys(char* argv[]) {
	int d1=0, d2=0;

	p = new polynomial*[2];
	if(!strcmp(argv[1],"-read")) {	/* if -read is given */

		if(!strcmp(argv[2],"-i")) {	/* read from file */
			ifstream myfile(argv[3]);
			if(myfile.is_open()){
				if(!strcmp(argv[4],"-d1"))
					d1 = atoi(argv[5]);
				if(!strcmp(argv[6],"-d2"))
					d2 = atoi(argv[7]);

				string line, str;

				for(int k=0; k<2; k++) {
					getline(myfile,line);
					cout << line << endl;

					str.clear();
					for(int i=0; i<line.size(); i++)
						if(line[i]!=' ')
							str.append(&(line[i]),1);	/* eliminate space */

					if(str[0]!='+')
						str.insert(0,"+");

					size_t pos, hold=0;
					while((pos = str.find_first_of("+-",hold)) != string::npos) {
						if(!isdigit(str[pos+1]))
							str.insert(pos+1,"1*");
						hold = pos+3;
					}

					str.erase(0,1);
					hold=0;
					while((pos = str.find("-",hold)) != string::npos) {
						str.insert(pos,"+");	/* insert "+" to separate terms */ 
						hold = pos+2;
					}


					int d = (k==0)?d1:d2;
					p[k] = new polynomial(d,str);
//					p[k]->print();
				}
				myfile.close();
			}
			else {
				cout << "Unable to open " << argv[3] << endl;
				exit(2);
			}
		}
		else if(!strcmp(argv[2],"-console")) {	/* read from console */
			if(!strcmp(argv[3],"-d1"))
				d1 = atoi(argv[4]);
			if(!strcmp(argv[5],"-d2"))
				d2 = atoi(argv[6]);

			string line, str;

			for(int k=0; k<2; k++) {
				getline(cin,line);
				cout << "Line is " << line << endl;

				str.clear();
				for(int i=0; i<line.size(); i++)
					if(line[i]!=' ')
						str.append(&(line[i]),1);	/* eliminate space */

				str.insert(0,"+");
				size_t pos, hold=0;
				while((pos = str.find_first_of("+-",hold)) != string::npos) {
					if(!isdigit(str[pos+1]))
						str.insert(pos+1,"1*");
					hold = pos+3;
				}

				str.erase(0,1);
				hold=0;
				while((pos = str.find("-",hold)) != string::npos) {
					str.insert(pos,"+");	/* insert "+" to separate terms */ 
					hold = pos+2;
				}

				int d = (k==0)?d1:d2;
				p[k] = new polynomial(d,str);
//				p[k]->print();
			}

		}

	}
	else if(!strcmp(argv[1],"-generate")) {
		if(!strcmp(argv[2],"-d1"))
			d1 = atoi(argv[3]);
		if(!strcmp(argv[4],"-d2"))
			d2 = atoi(argv[5]);

		srand(time(NULL));

		for(int i=0; i<2; i++) {
			int d = (i==0)?d1:d2;
			p[i] = new polynomial(d, "");     /* Empty string signifies make random polynomial*/
			cout << endl;
//			p[i]->print();
//			p[i]->printGenerate();
		}
	}
	else {
		cout << "Wrong arguments" << endl;
		exit(3);
	}

	/* Determine max(dx) and max(dy) for system */
	int d1x = p[0]->get_d('x');
	int d1y = p[0]->get_d('y');
	int d2x = p[1]->get_d('x');
	int d2y = p[1]->get_d('y');
//	cout << "d1x is " << d1x << " and d1y is " << d1y << " and d2x is " << d2x << " and d2y is " << d2y << endl;
	int maxdx = (d1x>d2x)?d1x:d2x;
	int maxdy = (d1y>d2y)?d1y:d2y;
//	this->hidden = (((d1y<d2y)?d1y:d2y) < ((d1x<d2x)?d1x:d2x) ) ? 'y':'x';
	this->hidden = (maxdy < maxdx)? 'y':'x';
	this->col = (this->hidden=='y')?maxdx+1:maxdy+1;
	this->depth = (this->hidden=='y')?maxdy+1:maxdx+1;
//	cout << "Maxdx is " << maxdx << " and maxdy is " << maxdy << " and hidden is " << hidden << " and col is " << col << " and depth is " << depth << endl;


	/* sys is a 3D matrix of size [2][max degree of non-hidden variable + 1][max degree of hidden varible + 1] */
	system_matrix = new double**[2];
	for(int i=0; i<2; i++) {
		system_matrix[i] = new double*[this->col];
		for(int j=0; j< this->col; j++) {
			system_matrix[i][j] = new double[this->depth];
		}
	}

	for(int i=0; i<2; i++)
		for(int j=0; j< this->col; j++)
			for(int k=0; k< this->depth; k++)
				system_matrix[i][j][k] = 0;

	/* Fill each vector with the correct vectors from the polyonym constant matrix */
	for(int i=0; i<2; i++)
		for(int j=0; j< this->col; j++)
			for(int k=0; k< this->depth; k++)
				p[i]->get_cons(system_matrix[i][j][k], j, k, this->hidden, this->depth);
}

sys::~sys() {

	for(int i=0; i<2; i++) {
		for(int j=0; j< this->col; j++) {
			delete [] system_matrix[i][j];
		}
		delete [] system_matrix[i];
	}
	delete [] system_matrix;

	delete p[0];
	delete p[1];
	delete [] p ;

}

void sys::print() {

	p[0]->print();
	p[1]->print();

	cout << endl;
	for(int i=0; i<2; i++) {
		for(int j=0; j<this->col; j++) {
			cout << "[";
			for(int k=0; k<this->depth-1; k++) {
				cout << system_matrix[i][j][k] << ",";
			}
			cout << system_matrix[i][j][this->depth-1]<< "]\t" ;
		}
		cout << endl;
	}
	cout << endl;
}

int sys::getDepth() {
		return this->depth;
}

/* Return a pointer to polynomial i */
polynomial* sys::get_pol(int i){
	if (i > -1 && i < 2)
		return this->p[i];
	else
		return NULL;
}

/* Return the degree of the non-hidden variable for polyonym pol */
int sys::get_d(int pol) {

	if (hidden == 'y'){
		return p[pol]->get_d('x');
	}
	else{
		return p[pol]->get_d('y');
	}
	return -1;
}

/* Write the pol[0] or pol[1] half of sys into the 2D matrix dest, starting at a specific row */
void sys::get_sys(double** dest, int pol, int skip){

	int columns;
	if (hidden == 'y')
		columns = p[pol]->get_d('x') + 1;
	else
		columns = p[pol]->get_d('y') + 1;

	int i, j;
	for (i = 0; i < columns; i++){
		for (j = 0; j < depth; j++){
			dest[i + skip][j] = system_matrix[pol][(columns - 1) - i][j];
		}
	}
	return;
	
}


