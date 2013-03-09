#ifndef SYMMETRY_SOLVER_H
#define SYMMETRY_SOLVER_H

#include <iostream>
#include <iomanip>
#include <string.h>
#include "symmetry_matrix.h"
#include "../unit_lower_triangular_matrix/ultriangular_matrix.h"
#include <ctime>

using std::cout;
using std::endl;

// #include <sys/time.h>

/*!	\brief Saves a permutation vector vec as a permutation matrix in matrix market (.mtx) format.
	\param vec the permutation vector.
	\param filename the filename the matrix will be saved under.
*/
template<class el_type>
bool save_perm(const std::vector<el_type>& vec, std::string filename)
{
	std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
	if(!out)
		return false;

	std::string header = "%%MatrixMarket matrix coordinate real general";;

	out << header << std::endl;
	out << vec.size() << " " << vec.size() << " " << vec.size() << "\n";

	for(int i = 0; i < (int) vec.size(); i++)
		out << i+1 << " " << vec[i]+1 << " " << 1 << "\n";
	
	out.close();
	return true;
}

/*! \brief Set of tools that facilitates conversion between different matrix formats. Also contains solver methods for matrices using a common interface.

	Currently, the only matrix type accepted is the lilc_matrix (as no other matrix type has been created yet).
*/
template<class el_type>
class solver
{
	public:
		symmetry_matrix<el_type> A;	///<The matrix to be factored.
		ultriangular_matrix<el_type> L;	///<The lower triangular factor of A.
		vector<int> perm;	///<A permutation vector containing all permutations on A.
		block_diag_matrix<el_type> D;	///<The diagonal factor of A.
		
		/*! \brief Loads the matrix A into solver.
			\param filename the filename of the matrix.
		*/
		void load(std::string filename)
		{
			//assert( A.load(filename) );
			clock_t start = clock(); double dif;
			A.load(filename);
			dif = clock() - start;
			printf("Load:\t%.3f seconds. \n", dif/CLOCKS_PER_SEC);
			printf("A is %d by %d with %d non-zeros.\n", A.n_rows(), A.n_cols(), A.nnz() );
		}
		
		/*! \brief Factors the matrix A into P' * S^(-1) * A * S^(-1) * P = LDL' in addition to printing some timing data to screen.
			
			More information about the parameters can be found in the documentation for the ildl() function.
			
			\param fill_factor a factor controling memory usage of factorization.
			\param tol a factor controling accuracy of factorization.
		*/
		void solve(double fill_factor, double tol)
		{
			perm.reserve(A.n_cols());
			// for (int i = 0; i < A.n_cols(); i++) {
				// perm.push_back(i);
			// }
			cout << std::fixed << std::setprecision(3);

			clock_t start = clock(); double dif, total = 0;
			
			A.equilibrate();
			
			dif = clock() - start; total += dif;
			printf("Equilibration:\t%.3f seconds.\n", dif/CLOCKS_PER_SEC);
			start = clock();
			
			A.rcm(perm);
			
			dif = clock() - start; total += dif;
			printf("RCM:\t\t%.3f seconds.\n", dif/CLOCKS_PER_SEC);
			
			start = clock();
			A.permute(perm);
			dif = clock() - start; total += dif;
			printf("Permutation:\t%.3f seconds.\n", dif/CLOCKS_PER_SEC);
			
			start = clock();
			A.ildl(L, D, perm, fill_factor, tol);
			dif = clock() - start; total += dif;
			
			printf("Factorization:\t%.3f seconds.\n", dif/CLOCKS_PER_SEC);	
			printf("Total time:\t%.3f seconds.\n", total/CLOCKS_PER_SEC);
			printf("L is %d by %d with %d non-zeros.\n", L.n_rows(), L.n_cols(), L.nnz() );

		}
		
		/*! \brief Save results of factorization (automatically saved into the output_matrices folder).
			
			The names of the output matrices follow the format out{}.mtx, where {} describes what the file contains (i.e. A, L, or D).
		*/
		void save()
		{
			cout << "Saving matrices..." << endl;
			A.save("output_matrices/outB.mtx", "output_matrices/outS.mtx");
			save_perm(perm, "output_matrices/outP.mtx");
			L.save("output_matrices/outL.mtx");
			D.save("output_matrices/outD.mtx");
			cout << "Save complete." << endl;
		}
		
		/*! \brief Prints the L and D factors to stdout.
		*/
		void display()
		{
			cout << A << endl;
			cout << L << endl;
			cout << D << endl;
			cout << "Permutation vector is: [" << perm << "]" << endl;
		}
};

#endif