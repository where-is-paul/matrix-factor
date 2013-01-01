//#include <iostream>
//#include <cassert>
//#include <string.h>
//#include "source/solver.h"
//#include "source/skew_solver.h"
//using namespace std;
//
//int main(int argc, char* argv[]) {
//
//	if (argc < 2 || argc > 6) {
//		std::cout << "Too many or too few arguments." << std::endl
//		<< "Program usage: ./ldl_driver [in.mtx] [lfil] [tol] [-save] [-display]" << std::endl;
//		return 0;
//	}
//
//	skew_solver<double> solv;
//	solv.load(argv[1]);
//	
//	double fill = 1.0, tol = 0.001;
//	if (argc > 2) fill = atof(argv[2]);
//	if (argc > 3) tol = atof(argv[3]);
//	
//	solv.solve(fill, tol);
//	
//	if (argc > 4) {
//		if (strcmp(argv[4], "-y") == 0) {
//			solv.save();
//		}
//		
//		if (argc > 5 && strcmp(argv[5], "-y") == 0) {
//			solv.display();
//		}
//		cout << "Done." << endl;
//	} else {
//		solv.save();
//	}
//
//	return 0;
//}