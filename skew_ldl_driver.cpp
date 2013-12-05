#include <iostream>
#include <cassert>
#include <string.h>
#include "source/skew_symmetric_matrix/skew_symmetry_solver.h"

using namespace std;

int main(int argc, char* argv[]) {

	if (argc < 2 || argc > 7) {
		std::cout << "Too many or too few arguments." << std::endl
		<< "Program usage: ./ldl_driver [in.mtx] [lfil] [tol] [-rcm/amd] [-save] [-display]" << std::endl;
		std::cout << "Sample usage: ./ldl_driver test_matrices/testmat1.mtx 1.0 0.001 -rcm -y -y" << endl;
		return 0;
	}

	skew_solver<double> solv;
	solv.load(argv[1]);
	
	double fill = 1.0, tol = 0.001;
	if (argc > 2) fill = atof(argv[2]);
	if (argc > 3) tol = atof(argv[3]);
	
	//default reordering scheme is AMD
	if (argc > 4) solv.set_reorder_scheme(argv[4]); 
	
	solv.solve_pp(fill, tol);
	
	if (argc > 5) {
		if (strcmp(argv[5], "-y") == 0) {
			solv.save();
		}
		
		if (argc > 6 && strcmp(argv[6], "-y") == 0) {
			solv.display();
			cout << endl;
		}
		cout << "Factorization Complete. All output written to /output_matrices directory." << endl;
	} else {
		solv.save();
	}

	return 0;
}