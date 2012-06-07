#include <iostream>
#include <string>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Crout.h>

using namespace Eigen;
using namespace std;

#define NDEBUG

int main(int argc, char* argv[]) {
    if (argc != 3) {
		cout << "Too many or too few arguments." << endl;
		cout << "Program usage: ./ldl_driver [in.mtx] [out.mtx]" << endl;
		return 0;
	}
	
	SparseMatrix<double> A;
	
	loadMarket(A, argv[1]);
	
	//SparseMatrix<double> L(A);
	VectorXd D(A.rows());
	
	CroutLDL<double>(A, D, false);
	
	saveMarket(A, argv[2]);
	
	//printf("The residual is %f.\n",((L * D.asDiagonal() * L.transpose()).triangularView<Lower>() - A).cwiseAbs().sum());
	
	return 0;
}
