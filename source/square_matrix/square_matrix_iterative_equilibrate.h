//-*- mode: c++ -*-
#ifndef SQUARE_MATRIX_ITERATIVE_EQUILIBRATE_H
#define SQUARE_MATRIX_ITERATIVE_EQUILIBRATE_H

using std::abs;
using std::min;
using std::max;

namespace iter_equil_helper {
	/*! \brief Computes the p-norm of column i of A. 
		\param A the input matrix.
		\param p the p-norm that we are using.
		\param i the column of the norm we wish to compute.
		\param temp storage space for extra information (the memory locations of elements with index < i) that this function computes.
		\return the p-norm of the i-th column of A.
	*/
	template<class el_type, class elt_it>
	inline double norm(square_matrix<el_type>* A, double p, int i, vector<elt_it>& temp) {
		typename square_matrix<el_type>::idx_it idx_iterator, i_it;
		typename square_matrix<el_type>::elt_it elt_iterator, e_it;
		double norm_val = 0;
		for (i_it = A->list[i].begin(); i_it != A->list[i].end(); i_it++) {
			A->coeffRef(i, *i_it, idx_iterator, elt_iterator);
			temp[*i_it] = elt_iterator;
			if (p < 0) {	//norm < 0 means INFINITY norm.
				norm_val = max(norm_val, abs(*elt_iterator));
			} else { //otherwise use the p norm.
				norm_val += pow(abs(*elt_iterator), p);
			}
		}

		for (e_it = A->m_x[i].begin(); e_it != A->m_x[i].end(); e_it++) {
			if (p < 0) {	//norm < 0 means INFINITY norm.
				norm_val = max(norm_val, abs(*e_it));
			} else { //otherwise use the p norm.
				norm_val += pow(abs(*e_it), p);
			}
		}
		
		if (p > 0) norm_val = pow(norm_val, 1.0/p);
		
		return norm_val;
	}
}

template<class el_type>
void square_matrix<el_type>::iterative_equilibrate(double norm, double tol, int max_iter) {
	// default norm is -1 (INF), tol is 1e-4
	//find termination points for loops with binary search later.
	int i, ncols = n_cols();
	idx_it idx_iterator;
	elt_it elt_iterator;
	vector<elt_it> temp(ncols);
	S.clear(); S.resize(ncols, 1);
	double norm_val, scale, min_norm = 10.0;
	// The algorithm proceeds in two stages: an inner iteration and an outer iteration.
	// 	- An inner iteration picks a particular row (say row k) and scales row k and col k
	// 	  by 1/sqrt(||row k||_inf).
	//	- An outer iteration applies n inner iterations, to row 1, 2, ... n.
	// Repeat outer iterations until convergence.
	
	// Note that after the first iteration, all values in the matrix will be <= 1.
	// To measure the convergence of this algorithm, we just check for when
	// abs(1 - min_i || row i ||_inf) goes under a specified tolerance.
	// (since the smallest row norm lower bounds all the other row norms).
	while (abs(1-min_norm) > tol && max_iter--) {
		//An outer inner iteration:
		for (i = 0; i < ncols; i++) {
			//assumes indices are ordered. since this procedure is run
			//before factorization pivots matrix, this is a fair assumption
			//for most matrix market matrices.
			
			//note that this is not the most efficient way of implementing
			//this equilibration procedure. we could instead keep an external vector
			//that maintains a partially calculated list of row/column norms.
			//however, equilibration is not the bottleneck in this program,
			//so this optimization is not so important (though I think it might
			//be able to speed up the equilibration by a factor of 2)
			
			//determine the row norm of column i:
			norm_val = iter_equil_helper::norm(this, norm, i, temp);
			
			//scale row i and column i by 1/sqrt(norm)
			scale = 1.0/sqrt(norm_val);
			if (scale > eps) {
				S[i] *= scale;
				for (idx_it it = list[i].begin(); it != list[i].end(); it++)
					//can use bin. search on coeff since no reordering is done yet.
					*temp[*it] *= scale;

				if (!m_idx[i].empty() && (m_idx[i][0] == i))
					m_x[i][0] *= scale;
				for (elt_it it = m_x[i].begin(); it != m_x[i].end(); it++)
					*it *= scale;
			}
		}
		
		// find minimum norm of over all columns
		min_norm = 10.0;
		for (int i = 0; i < ncols; i++) {
			norm_val = iter_equil_helper::norm(this, norm, i, temp);
			min_norm = min(min_norm, norm_val);
		}
		
		//printf("%f\n", min_norm);
		//fflush(stdout);
	}

}

#endif