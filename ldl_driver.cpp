#include <iostream>
#include <cassert>
#include "csc_matrix.h"

using namespace std;

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cout << "Too many or too few arguments." << std::endl
				  << "Program usage: ./ldl_driver [in.mtx]" << std::endl;
		return 0;
	}
	
	csc_matrix<int, double> A, L;
	vector<double> D;
	
	assert( A.load(argv[1]) );
	//cout << A << endl;
	
	A.ildl(L, D, 30, 0.00);
	cout << L << endl;
	cout << D << endl;	
	
	return 0;
}
