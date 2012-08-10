#include <iostream>
#include <sys/time.h>
#include <string.h>
#include "lilc_matrix.h"

template<class el_type>
bool save_perm(const std::vector<el_type>& vec, std::string filename) {
	std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
	if(!out)
	return false;

	out.flags(std::ios_base::scientific);
	out.precision(12);
	std::string header = "%%MatrixMarket matrix coordinate real general";; 

	out << header << std::endl; 
	out << vec.size() << " " << vec.size() << " " << vec.size() << "\n";

	for(int i = 0; i < (int) vec.size(); i++)
	out << i+1 << " " << i+1 << " " << vec[i]+1 << "\n";
	
	out.close();
	return true;
}

/*! \brief Set of tools that facilitates conversion between different matrix formats. Also contains solver methods for matrices using a common interface.

	Currently, the only matrix type accepted is the lilc_matrix (as no other matrix type has been created yet).
*/
template<class el_type, class mat_type = lilc_matrix<el_type> >
class solver
{
	public:
		mat_type A, L;
		vector<int> perm;
		block_diag_matrix<el_type> D;
		
		void load(std::string filename) {
			assert( A.load(filename) );
			printf("A is %d by %d with %d non-zeros.\n", A.n_rows(), A.n_cols(), A.nnz() );
		}
		
		void solve(double fill_factor, double tol) {
			perm.reserve(A.n_cols());
			
			A.sym_equil();
			A.sym_rcm(perm);
			A.sym_perm(perm);
			
			struct timeval tim;  
			gettimeofday(&tim, NULL);  
			double t1=tim.tv_sec+(tim.tv_usec/1e6);  
			
			A.ildl(L, D, perm, fill_factor, tol);
			
			gettimeofday(&tim, NULL);  
			double t2=tim.tv_sec+(tim.tv_usec/1e6);  
			
			printf("The factorization took %.6lf seconds.\n", t2-t1);
			printf("L is %d by %d with %d non-zeros.\n", L.n_rows(), L.n_cols(), L.nnz() ); 
		}
		
		void save() {
			cout << "Saving matrices..." << endl;
			A.save("output_matrices/outA.mtx", true);
			A.S.save("output_matrices/outS.mtx");
			save_perm(perm, "output_matrices/outPerm.mtx");
			L.save("output_matrices/outL.mtx");
			D.save("output_matrices/outD.mtx");
		}
		
		void display() {
			cout << L << endl;
			cout << D << endl;
			cout << perm << endl;
		}
};