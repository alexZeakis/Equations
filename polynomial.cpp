#include <iostream>
#include <stdlib.h>
#include "polynomial.h"

#define RANGE 100


polynomial::polynomial(int d, string pol) {
	this->d = d;

	cons = new int*[d+1];
	for(int i=0; i<d+1; i++)
		cons[i] = new int[d+1];

	for(int i=0; i<d+1; i++)
		for(int j=0; j<d+1; j++)
			cons[i][j]= 0;	/* initializing */

	if (pol.size() > 0) {                           /* Create polyonym from user input */
		size_t pos =0;
		while(pol.size() > 0) {
			string token;

			if((pos = pol.find("+")) != string::npos) {
				token = pol.substr(0,pos);
			}
			else {	/* last term in string */
				token = pol;
				pos = pol.size();
			}

			string c,var1,var2;
			int con, degx=0, degy=0;
			size_t dpos;
			if((dpos = token.find("*")) != string::npos) {	/* Check if term that includes variables or constant */
				c = token.substr(0,dpos);

				con = atoi(c.c_str());
				token.erase(0,dpos+1);

				if ((dpos = token.find("*")) != string::npos) {	/* Check if term with two variables or one */
					var1 = token.substr(0,dpos);
					var2 = token.substr(dpos+1,token.size());
					if(var1[0]=='x') {
							var1.erase(0,2);	
						var2.erase(0,2);	
						degx = atoi(var1.c_str());           /* Save the exponents for x and y*/
						degy = atoi(var2.c_str());
					}
					else { /* var1[0]=='y' */
						var1.erase(0,2);	
						var2.erase(0,2);	
						degy = atoi(var1.c_str());
						degx = atoi(var2.c_str());
					}
					if(!degx) degx=1;
					if(!degy) degy=1;
				}
				else { /* term with single variable */
					if(token[0] == 'x') {
						token.erase(0,2);
						degx = atoi(token.c_str());
						if(!degx) degx=1;
					}
					else { /* token[0]=='y' */
						token.erase(0,2);
						degy = atoi(token.c_str());
						if(!degy) degy=1;
					}
				}
			}

			else {
				con = atoi(token.c_str());
			}
			pol.erase(0, pos+1);
			cons[degx][degy] = con;
		}
	}
	else {		/* Create polyonym randomly by generating random constants for the upper left half of the 
				2D cons matrix */ 
	
		int non_zero = 0;         /* flag for non-zero values on the anti-digonal of the matrix */
		for (int i = 0; i < d + 1; i++){
			for (int j = 0; j < (d + 1) - i; j++){
				int c = (rand() % (RANGE + 1)) - RANGE / 2;
				if (!non_zero && j == d - i){
					if (c != 0)
						non_zero = 1;
					else if (i == d){  /* if the rest of the anti-diagonal is zero, cons[d][0] must be non-zero*/
						while (c == 0)
							c = (rand() % (RANGE + 1)) - RANGE / 2;
					}
				}
				cons[i][j] = c;
			}
		}

	}
}


polynomial::~polynomial() {

	for (int i = 0; i < d + 1; i++)
		delete [] cons[i];
	delete [] cons;
}


void polynomial::print() {
	cout << endl;
	for(int i=0; i<d+1; i++) {
		for(int j=0; j<d+1; j++)
			cout << cons[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
}

/* Return dx for get_d('x'), dy for get_d('y') */
int polynomial::get_d(char var){

	int non_zero = 0;
	int i, j;

	if (var != 'x' && var != 'y')
		return -1;
	
	/*Look for a non-zero value in the cons matrix starting from bottom row (for x) or last column (for y)*/
	for (i = d; i >= 0; i--){
		for (j = 0; j < d + 1; j++){
			if ((var == 'x' && cons[i][j] != 0) || (var == 'y' && cons[j][i] != 0)){
				non_zero = 1;
				break;
			}
		}
		if (non_zero)
			break;
	}
	return i;

}

/* Write the vector of constants for a specific degree of x or a specific degree of y into vector dest */
void polynomial::get_cons(int& dest, int j, int k, char var, int depth){

	if(j< d+1) {
		if (var == 'y'){
			dest = cons[j][k];
		}
		else if (var == 'x'){
			dest = cons[k][j];
		}
	}

}

