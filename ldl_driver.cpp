#include <iostream>
#include <cassert>
#include <string.h>
#include <sys/time.h>
#include "csc_matrix.h"

using namespace std;

/*! \mainpage Main Page
 *
 * \section intro_sec Introduction
 *
 * This page will contain some information on how to use the package in the future. For now it is empty.
 *
 */
 
int main(int argc, char* argv[]) {
	if (argc < 4) {
		std::cout << "Too many or too few arguments." << std::endl
				  << "Program usage: ./ldl_driver [lfil] [tol] [in.mtx] [L_out.mtx] [D_out.mtx] [-yn]" << std::endl;
		return 0;
	}
	
	csc_matrix<int, double> A, L;
	vector<double> D;
	
	assert( A.load(argv[3]) );
	//cout << A << endl;
	
	struct timeval tim;  
    gettimeofday(&tim, NULL);  
    double t1=tim.tv_sec+(tim.tv_usec/1e6);  
    
	A.ildl(L, D, atof(argv[1]), atof(argv[2]));
	
	gettimeofday(&tim, NULL);  
    double t2=tim.tv_sec+(tim.tv_usec/1e6);  
    printf("The factorization took %.6lf seconds.\n", t2-t1);
	printf("L has %d non-zeros.\n", L.nnz()); 
	
	if (argc > 4) {
		L.save(argv[4]);
		save(D, argv[5]);
		if (argc > 6 && strcmp(argv[6], "-y") == 0) {
			cout << L << endl;
			cout << D << endl;
		}
	}
	
	return 0;
}
