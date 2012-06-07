#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

// this is just a test file of eigen's capabilities
int main() {
	
	SparseMatrix<double> A;
	
	loadMarket(A, "aug3dcqp.mtx");
    IncompleteLUT<double> ILUT(A, 0.01, 30);
	//ILUT.compute(A);
          
    return 0;
}
