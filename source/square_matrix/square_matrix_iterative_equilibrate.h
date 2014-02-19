//-*- mode: c++ -*-
#ifndef SQUARE_MATRIX_EQUILIBRATE_H
#define SQUARE_MATRIX_EQUILIBRATE_H

using std::abs;

namespace iter_equi_helper {
	template<class el_type, class elt_it>
	double norm(square_matrix<el_type>& A, int i, vector<elt_it>& temp) {
		square_matrix<el_type>::idx_it idx_iterator;
		square_matrix<el_type>::elt_it elt_iterator;
		double norm_val = 0;
		for (square_matrix<el_type>::idx_it it = list[i].begin(); it != list[i].end(); it++) {
			coeffRef(i, *it, idx_iterator, elt_iterator);
			temp[*it] = elt_iterator;
			if (norm < 0) {	//norm < 0 means INFINITY norm.
				norm_val = max(abs(*elt_iterator));
			} else { //otherwise use the p norm.
				norm_val += pow(abs(*elt_iterator), norm);
			}
		}

		for (square_matrx<el_type>::elt_it it = m_x[i].begin(); it != m_x[i].end(); it++) {
			if (norm < 0) {	//norm < 0 means INFINITY norm.
				norm_val = max(abs(*it));
			} else { //otherwise use the p norm.
				norm_val += pow(abs(*it), norm);
			}
		}
		
		if (norm > 0) norm_val = pow(norm_val, 1.0/norm);
		
		return norm;
	}
}

template<class el_type>
void square_matrix<el_type>::iterative_equilibrate(double norm = -1, double tol = 1e-4)
{
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
	while (abs(1-min_norm) > tol) {
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
			norm_val = iter_equil_helper::norm(this, i, temp);
			
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
		for (int i = 0; i < ncols; i++) {
			norm_val = iter_equil_helper::norm(this, i, temp);
			min_norm = min(min_norm, norm_val);
		}
	}

}

#endif