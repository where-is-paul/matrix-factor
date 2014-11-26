#ifndef SKEW_SYMMETRY_MATRIX_ILDLRP_H
#define SKEW_SYMMETRY_MATRIX_ILDLRP_H

using std::abs;
using std::cout;
using std::min;
using std::endl;

template <class el_type>
void skew_symmetry_matrix<el_type>::ildlrp(ultriangular_matrix<el_type>& L, skew_block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol, const double& pp_tol)
{
	//----------------- initialize temporary variables --------------------//
	const int ncols = n_cols(); // number of cols in A.
	const double stat_piv = 1e-6;

	int lfil;
	if (fill_factor >= 1e4) lfil = ncols; // just incase users decide to enter a giant fill factor for fun...
	else lfil = (int) (fill_factor*nnz()/ncols); // roughly a factor of 2 since only lower tri. of A is stored
	
	el_type wi, wr;
	int i, j, k, r, newr;

	vector<bool> in_set(ncols, false); // bitset used for unsorted merges
	square_matrix_swap_struct<el_type> s;	// struct containing temp vars used in pivoting.

	elt_vector_type col_i(ncols, 0), col_r(ncols, 0); // work vector for the current columns
	idx_vector_type col_i_nnzs, col_r_nnzs;  // non-zeros on current col.
	col_i_nnzs.reserve(ncols); // reserves space for worse case (entire col is non-zero)
	col_r_nnzs.reserve(ncols);

	//--------------- allocate memory for L and D ------------------//
	L.resize(ncols, ncols, min(2*lfil, 32)); //allocate a vector of size n for Llist as well
	D.resize(ncols );

	//------------------- main loop: factoring begins -------------------------//
	for (k = 0; k < ncols; k=k+2)
	{
		cout << k << endl;
		i = k;
		// assign nonzeros indices and values of A(k:n, i) to col_i_nnzs
		col_i_nnzs.assign(m_idx[k].begin(), m_idx[k].end());
		for (j = 0; j < (int) col_i_nnzs.size(); j++)
			col_i[col_i_nnzs[j]] = m_x[k][j];
		// do delayed updates on current column
		// A(k:n, r) += L(k:n, i+1)*D(i, i+1)*L(r, i) + L(k:n, i)*D(i+1, i)*L(r, i+1) (i = 0, 2, ..., k-2)
		update(i, col_i, col_i_nnzs, L, D, in_set);
		wi = max(col_i, col_i_nnzs, r);

		while (true)
		{
			// assign nonzeros indices and values of A(k:n, r) to col_r_nnzs
			for (idx_it it = col_r_nnzs.begin(); it != col_r_nnzs.end(); it++)
				col_r[*it] = 0;
			col_r_nnzs.clear();
			for (j = list_first[r]; j < (int) list[r].size(); j++)
			{
				col_r_nnzs.push_back(list[r][j]);
				col_r[list[r][j]] = -coeff(r, list[r][j]);
			}
			col_r_nnzs.insert(col_r_nnzs.end(), m_idx[r].begin(), m_idx[r].end());
			for (j = 0; j < (int) m_idx[r].size(); j++)
				col_r[m_idx[r][j]] = m_x[r][j];
			update(r, col_r, col_r_nnzs, L, D, in_set);

			//find maximum element in work and store its index in r.
			wr = max(col_r, col_r_nnzs, newr);
			
			if (abs(wi - wr) < eps)
			{
				// swap rows and columns k and i, k+1 and r
				if (k != i)
				{
					pivot(s, L, k, i);
					std::swap(perm[k], perm[i]);
					std::swap(col_i[k], col_i[i]);
					std::swap(col_r[k], col_r[i]);
					safe_swap(col_i_nnzs, k, i);
					safe_swap(col_r_nnzs, k, i);
				}

				advance_list(k);
				L.advance_column(k);

				if (k+1 != r)
				{
					pivot(s, L, k+1, r);
					std::swap(perm[k+1], perm[r]);
					std::swap(col_i[k+1], col_i[r]);
					std::swap(col_r[k+1], col_r[r]);
					safe_swap(col_i_nnzs, k+1, r);
					safe_swap(col_r_nnzs, k+1, r);
				}
				break;
			}
			else
			{
				i = r;
				wi = wr;
				r = newr;
				col_i.swap(col_r);
				col_i_nnzs.swap(col_r_nnzs);
			}
		}
		//--------------end pivoting--------------//

		calculate(k, L, D, col_i, col_r, col_i_nnzs, col_r_nnzs, tol, lfil, stat_piv);
	}

	L.nnz_count += ncols; // diagonal 1s are not stored in L, but L.nnz_count should include them
}

#endif
