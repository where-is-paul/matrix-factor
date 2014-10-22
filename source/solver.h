#ifndef _SOLVER_H
#define _SOLVER_H

#include <iostream>
#include <string.h>
#include "lilc_matrix.h"
#include <ctime>
#include <iomanip>

/*!	\brief Saves a permutation vector vec as a permutation matrix in matrix market (.mtx) format.
  \param vec the permutation vector.
  \param filename the filename the matrix will be saved under.
 */
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

	for(int i = 0; i < (int) vec.size(); i++) {
		out << i+1 << " " << 1 << " " << vec[i]+1 << "\n";
	}
	out.close();
	return true;
}

/*! \brief Set of tools that facilitates conversion between different matrix formats. Also contains solver methods for matrices using a common interface.

  Currently, the only matrix type accepted is the lilc_matrix (as no other matrix type has been created yet).
 */
template<class el_type, class mat_type = lilc_matrix<el_type> >
class solver {
	public:
		mat_type A;	///<The matrix to be factored.
		mat_type L;	///<The lower triangular factor of A.
		vector<int> perm;	///<A permutation vector containing all permutations on A.
		block_diag_matrix<el_type> D;	///<The diagonal factor of A.
		int reorder_scheme; ///<Set to to 0 for AMD, 1 for RCM, 2 for no reordering.
		bool equil; ///<Set to true for max-norm equilibriation.
		bool has_rhs; ///<Set to true if we have a right hand side that we expect to solve.
		vector<el_type> rhs; ///<The right hand side we'll solve for.
		vector<el_type> sol_vec; ///<The solution vector.
		
		/*! \brief Solver constructor, initializes default reordering scheme.
		*/
		solver() {
		reorder_scheme = 0;
		equil = true;
	}

	/*! \brief Loads the matrix A into solver.
	  \param filename the filename of the matrix.
	 */
	void load(std::string filename) {
		bool result = A.load(filename);
		assert(result);
		printf("A is %d by %d with %d non-zeros.\n", A.n_rows(), A.n_cols(), A.nnz() );
	}

	/*! \brief Sets the reordering scheme for the solver.
	 */
	void set_reorder_scheme(const char* ordering) {
		if (strcmp(ordering, "rcm") == 0) {
			reorder_scheme = 1;
		} else if (strcmp(ordering, "amd") == 0) {
>>>>>>> 75c2fd562de7c8cf72002762a357cedab10988fe
			reorder_scheme = 0;
		} else if (strcmp(ordering, "none") == 0) {
			reorder_scheme = 2;
		}
		
		/*! \brief Loads a right hand side b into the solver.
			\param b a vector of the right hand side.
		*/
		void set_rhs(vector<el_type> b) {
			rhs = b;
			has_rhs = true;
			printf("Right hand side has %d entries.\n", rhs.size() );
		}
		
		/*! \brief Sets the reordering scheme for the solver.
		*/
		void set_reorder_scheme(const char* ordering) {
			if (strcmp(ordering, "rcm") == 0) {
					reorder_scheme = 1;
			} else if (strcmp(ordering, "amd") == 0) {
					reorder_scheme = 0;
			} else if (strcmp(ordering, "none") == 0) {
					reorder_scheme = 2;
			}
		}

	/*! \brief Factors the matrix A into P' * S^(-1) * A * S^(-1) * P = LDL' in addition to printing some timing data to screen.

	  More information about the parameters can be found in the documentation for the ildl() function.

	  \param fill_factor a factor controling memory usage of factorization.
	  \param tol a factor controling accuracy of factorization.
	 */
	void solve(double fill_factor, double tol, double pp_tol) {
		perm.reserve(A.n_cols());
		cout << std::fixed << std::setprecision(3);
		//gettimeofday(&tim, NULL);  
		//double t0=tim.tv_sec+(tim.tv_usec/1e6);
		clock_t start = clock(); double dif, total = 0;

		if (equil == 1) {
			A.sym_equil();
			dif = clock() - start; total += dif; 
			printf("Equilibriation:\t%.3f seconds.\n", dif/CLOCKS_PER_SEC);
		}
<<<<<<< HEAD
		
		/*! \brief Factors the matrix A into P' * S^(-1) * A * S^(-1) * P = LDL' in addition to printing some timing data to screen.
			
			More information about the parameters can be found in the documentation for the ildl() function.
			
			\param fill_factor a factor controling memory usage of factorization.
			\param tol a factor controling accuracy of factorization.
			\param pp_tol a factor controling the aggresiveness of Bunch-Kaufman pivoting.
			\param max_iter the maximum number of iterations for minres (ignored if no right hand side).
		*/
		void solve(double fill_factor, double tol, double pp_tol, int max_iter = -1) {
			perm.reserve(A.n_cols());
			cout << std::fixed << std::setprecision(3);
			//gettimeofday(&tim, NULL);  
			//double t0=tim.tv_sec+(tim.tv_usec/1e6);
			clock_t start = clock(); double dif, total = 0;

			if (equil == 1) {
				A.sym_equil();
				dif = clock() - start; total += dif; 
				printf("Equilibriation:\t%.3f seconds.\n", dif/CLOCKS_PER_SEC);
			}
=======
>>>>>>> 75c2fd562de7c8cf72002762a357cedab10988fe

		if (reorder_scheme != 2) {
			start = clock();
			std::string perm_name;
			switch (reorder_scheme) {
				case 0:
					A.sym_amd(perm);
					perm_name = "AMD";
					break;
				case 1:
					A.sym_rcm(perm);
					perm_name = "RCM";
					break;
			}

			dif = clock() - start; total += dif;
			printf("%s:\t\t%.3f seconds.\n", perm_name.c_str(), dif/CLOCKS_PER_SEC);

			start = clock();
			A.sym_perm(perm);
			dif = clock() - start; total += dif;
<<<<<<< HEAD
			
			printf("Factorization:\t%.3f seconds.\n", dif/CLOCKS_PER_SEC);
			printf("Total time:\t%.3f seconds.\n", total/CLOCKS_PER_SEC);
			printf("L is %d by %d with %d non-zeros.\n", L.n_rows(), L.n_cols(), L.nnz() );
			fflush(stdout);
			
			if (has_rhs) {
				assert(max_iter >= 0);
				minres(max_iter);
			}
		}
		
		/*! \brief Applies minres on A, preconditioning with factors L and D..
			
			\param max_iter the maximum number of minres iterations.
			\param stop_tol the stopping tolerance of minres. i.e. we stop as soon as the residual goes below stop_tol.
			\param shift shifts A by shift*(identity matrix) to make it more positive definite. This sometimes helps.
		*/
		void minres(int max_iter = 1000, double stop_tol = 1e-6, double shift = 0.0);
		
		/*! \brief Save results of factorization (automatically saved into the output_matrices folder).
			
			The names of the output matrices follow the format out{}.mtx, where {} describes what the file contains (i.e. A, L, or D).
		*/
		void save() {
			cout << "Saving matrices..." << endl;
			A.save("output_matrices/outB.mtx", true);
			A.S.save("output_matrices/outS.mtx");
			save_perm(perm, "output_matrices/outP.mtx");
			L.save("output_matrices/outL.mtx", false);
			D.save("output_matrices/outD.mtx");
			cout << "Save complete." << endl;
		}
		
		/*! \brief Prints the L and D factors to stdout.
		*/
		void display() {
			cout << L << endl;
			cout << D << endl;
			cout << perm << endl;
=======
			printf("Permutation:\t%.3f seconds.\n", dif/CLOCKS_PER_SEC);
		} else {
			// no permutation specified, store identity permutation instead.
			for (int i = 0; i < A.n_cols(); i++) {
				perm.push_back(i);
			}
>>>>>>> 75c2fd562de7c8cf72002762a357cedab10988fe
		}

		start = clock();
		A.ildl(L, D, perm, fill_factor, tol, pp_tol);
		dif = clock() - start; total += dif;

		printf("Factorization:\t%.3f seconds.\n", dif/CLOCKS_PER_SEC);
		printf("Total time:\t%.3f seconds.\n", total/CLOCKS_PER_SEC);
		printf("L is %d by %d with %d non-zeros.\n", L.n_rows(), L.n_cols(), L.nnz() );
		fflush(stdout);
	}

	/*! \brief Save results of factorization (automatically saved into the output_matrices folder).

	  The names of the output matrices follow the format out{}.mtx, where {} describes what the file contains (i.e. A, L, or D).
	 */
	void save() {
		cout << "Saving matrices..." << endl;
		A.save("output_matrices/outB.mtx", true);
		A.S.save("output_matrices/outS.mtx");
		save_perm(perm, "output_matrices/outP.mtx");
		L.save("output_matrices/outL.mtx", false);
		D.save("output_matrices/outD.mtx");
		cout << "Save complete." << endl;
	}

	/*! \brief Prints the L and D factors to stdout.
	 */
	void display() {
		cout << L << endl;
		cout << D << endl;
		cout << perm << endl;
	}
};

<<<<<<< HEAD
#include "solver_minres.h"

#endif
=======
#endif
>>>>>>> 75c2fd562de7c8cf72002762a357cedab10988fe
