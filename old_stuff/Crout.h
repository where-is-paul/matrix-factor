#include <Eigen/Core>
#include <iostream>

using namespace Eigen;
using namespace std;

//TODO: Add symmetry check and standard stuff to prevent foolishness.

template <typename DerivedA>
void CroutLDL(SparseMatrix<DerivedA> &L, 
              VectorXd &D, 
              bool pivot) {
              
    //Do some stuff...
	cout << "Beginning factorization..." << endl;
    int n = L.rows();
	
    if (pivot) {
        //Do some pivoting...
		//cout << "L:\n" << L << endl;
		//cout << L.subcols(2, 5) << endl;
    }	

	for (int k = 0; k < n; k++) {
		D(k) = L.coeff(k,k);
	}
	
	SparseVector<DerivedA> z(L.rows());
	for (int k = 0; k < n; k++) {
		//cout << k << endl;
		for (int i = 0; i < k; i++) {
			if (L.coeff(k,i) != 0) { //TODO: change this to using machine EPS later.
				z = L.coeff(k,i) * D(i) * L.col(i).eval();
				for (typename SparseVector<DerivedA>::InnerIterator it(z); it; ++it) {
					if (it.index() > k) {
						L.coeffRef(it.index(), k) -= it.value();
					}
				}
			}
		}
		
		L.col(k) /= D(k); L.coeffRef(k,k) = 1;
		for (typename SparseMatrix<DerivedA>::InnerIterator it(L, k); it; ++it) {
			if (it.row() > k)
				D(it.row()) -= it.value() * D(k) * it.value(); 
		}
	}
	
	cout << "Factorization complete" << endl;
}
