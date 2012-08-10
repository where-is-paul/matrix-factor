#include <iostream>
#include <cassert>
#include <string.h>
#include "source/solver.h"

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

	solver<double> solv;
	solv.load(argv[3]);
	solv.solve(atof(argv[1]), atof(argv[2]));
	
	if (argc > 4) {
		if (strcmp(argv[4], "-y") == 0) {
			solv.save();
		}
		
		if (argc > 5 && strcmp(argv[5], "-y") == 0) {
			solv.display();
		}
		cout << "Done." << endl;
	}

	return 0;
}
