#include <iostream>
#include <cassert>
#include <string.h>
#include <sys/time.h>
#include "source/lilc_matrix.h"

using namespace std;

/*! \mainpage Main Page
*
* \section intro_sec Introduction
*
* This page will contain some information on how to use the package in the future. For now it is empty.
*
*/

int main(int argc, char* argv[]) {

	if (argc < 4 || argc > 7) {
		std::cout << "Too many or too few arguments." << std::endl
		<< "Program usage: ./ldl_driver [lfil] [tol] [in.mtx] [-save] [-yn]" << std::endl;
		return 0;
	}

	
	lilc_matrix<double> A, L;
	vector<int> perm;
	block_diag_matrix<double> D;
	
	assert( A.load(argv[3]));
	printf("A is %d by %d with %d non-zeros.\n", A.n_rows(), A.n_cols(), A.nnz() );

	struct timeval tim;  
	gettimeofday(&tim, NULL);  
	double t1=tim.tv_sec+(tim.tv_usec/1e6);  
	
	perm.resize(A.n_cols());
	for (int i = 0; i < A.n_cols(); i++) perm[i] = i;
	
	A.ildl(L, D, perm, atof(argv[1]), atof(argv[2]));
	
	gettimeofday(&tim, NULL);  
	double t2=tim.tv_sec+(tim.tv_usec/1e6);  
	printf("The factorization took %.6lf seconds.\n", t2-t1);
	printf("L is %d by %d with %d non-zeros.\n", L.n_rows(), L.n_cols(), L.nnz() ); 
	
	if (argc > 4) {
		if (strcmp(argv[4], "-y") == 0) {
			cout << "Saving matrices..." << endl;
			L.save("output_matrices/outL.mtx");
			D.save("output_matrices/outD.mtx");
			save(perm, "output_matrices/outPerm.mtx");
		}
		
		if (argc > 5 && strcmp(argv[5], "-y") == 0) {
			cout << L << endl;
			cout << D << endl;
			cout << perm << endl;
		}
		cout << "Done." << endl;
	}

	return 0;
}
