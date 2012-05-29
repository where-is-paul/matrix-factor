#include <Eigen/Sparse>
#include <iostream>

using namespace Eigen;
using namespace std;

template <typename Derived>
void CroutLDL(SparseMatrix<Derived> &A, 
              SparseMatrix<Derived> &L, 
              DiagonalMatrix<Derived> &D, 
              bool pivot) {
              
    //Do some stuff...
    
    if (pivot) {
        //Do some pivoting...
    }
}
