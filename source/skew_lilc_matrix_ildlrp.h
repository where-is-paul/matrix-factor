#ifndef _SKEW_LILC_MATRIX_ILDLRP_H_
#define _SKEW_LILC_MATRIX_ILDLRP_H_

template <class el_type>
void skew_lilc_matrix<el_type>::ildlrp(lilc_matrix<el_type>& L, skew_block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol)
{
	//----------------- initialize temporary variables --------------------//
	const int ncols = n_cols(); // number of cols in A.
	int lfil = (int)(2 * fill_factor * nnz() / ncols); // roughly a factor of 2 since only lower tri. of A is stored
	el_type wi, wr, d;
	int count = 0; // the total number of nonzeros stored in L.
	int i, j, k, r, newr;
	vector<bool> in_set(ncols, false); // bitset used for unsorted merges
	skew_swap_struct<el_type> s;	// struct containing temp vars used in pivoting.
	elt_vector_type col_i(ncols, 0), col_r(ncols, 0); // work vector for the current columns
	idx_vector_type col_i_nnzs, col_r_nnzs;  // non-zeros on current col.
	col_i_nnzs.reserve(ncols); // reserves space for worse case (entire col is non-zero)
	col_r_nnzs.reserve(ncols);

	//--------------- allocate memory for L and D ------------------//
	L.resize(ncols, ncols);
	D.resize(ncols );

	//------------------- main loop: factoring begins -------------------------//
	for (k = 0; k < ncols; k=k+2)
	{
		i = k;
		// assign nonzeros indices and values of A(k:n, i) to col_i_nnzs
		col_i_nnzs.assign(m_idx[k].begin(), m_idx[k].end());
		for (j = 0; j < (int) col_i_nnzs.size(); j++)
			col_i[col_i_nnzs[j]] = m_x[k][j];
		// do delayed updates on current column
		// A(k:n, r) += L(k:n, i+1)*D(i, i+1)*L(r, i) + L(k:n, i)*D(i+1, i)*L(r, i+1) (i = 0, 2, ..., k-2)
		skew_update(i, col_i, col_i_nnzs, L, D, in_set);
		wi = max(col_i, col_i_nnzs, r);

		while (true)
		{
			// assign nonzeros indices and values of A(k:n, r) to col_r_nnzs
			for (auto it = col_r_nnzs.begin(); it != col_r_nnzs.end(); it++)
				col_r[*it] = 0;
			col_r_nnzs.clear();
			for (j = first[r]; j < (int) list[r].size(); j++)
			{
				col_r_nnzs.push_back(list[r][j]);
				col_r[list[r][j]] = -coeff(r, list[r][j]);
			}
			col_r_nnzs.insert(col_r_nnzs.end(), m_idx[r].begin(), m_idx[r].end());
			for (j = 0; j < (int) m_idx[r].size(); j++)
				col_r[m_idx[r][j]] = m_x[r][j];
			skew_update(r, col_r, col_r_nnzs, L, D, in_set);

			//find maximum element in work and store its index in r.
			wr = max(col_r, col_r_nnzs, newr);
			
			if (abs(wi - wr) < eps)
			{
				// swap rows and columns k and i, k+1 and r
				if (k != i)
				{
					pivot(s, in_set, L, k, i);
					std::swap(perm[k], perm[i]);
					std::swap(col_i[k], col_i[i]);
					std::swap(col_r[k], col_r[i]);
					safe_swap(col_i_nnzs, k, i);
					safe_swap(col_r_nnzs, k, i);
				}

				advance_list(k);
				L.advance_first(k);

				if (k+1 != r)
				{
					pivot(s, in_set, L, k+1, r);
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

		D[k] = d = col_i[k+1];
		col_i[k] = col_i[k+1] = col_r[k] = col_r[k+1] = 0;

		// erase diagonal, sub-diagonal, and super-diagonal elements from non-zero indices (to exclude it from being dropped)
		col_i_nnzs.erase(std::remove_if(col_i_nnzs.begin(), col_i_nnzs.end(), isKorKp1(k)), col_i_nnzs.end());
		col_r_nnzs.erase(std::remove_if(col_r_nnzs.begin(), col_r_nnzs.end(), isKorKp1(k)), col_r_nnzs.end());
		// performs the dual dropping procedure.
		drop_tol(col_i, col_i_nnzs, lfil, tol);
		drop_tol(col_r, col_r_nnzs, lfil, tol);
		
		// resize kth column of L to proper size.
		L.m_idx[k].resize(col_r_nnzs.size()+1);
		L.m_x[k].resize(col_r_nnzs.size()+1);
		L.m_idx[k][0] = k;
		L.m_x[k][0] = 1;
		
		if (abs(D[k]) < eps) D[k] = 1e-6; //statically pivot
		// fill L(k+2:n, k)
		i = 1;
		for (auto it = col_r_nnzs.begin(); it != col_r_nnzs.end(); it++)
		{
			if (abs(col_r[*it]) > eps)
			{
				L.m_idx[k][i] = *it; //col k nonzero indices of L are stored
				L.m_x[k][i] = -col_r[*it] / d; //col k nonzero values of L are stored
				L.list[*it].push_back(k); //update Llist
				count++;
				i++;
			}
		}
		// resize again, some 0s are dropped
		L.m_idx[k].resize(i);
		L.m_x[k].resize(i);


		//resize (k+1)th column of L to proper size.
		L.m_idx[k+1].resize(col_i_nnzs.size()+1);
		L.m_x[k+1].resize(col_i_nnzs.size()+1);
		L.m_idx[k+1][0] = k+1;
		L.m_x[k+1][0] = 1;
		
		// fill L(k+2:n, k+1)
		i = 1;
		for (auto it = col_i_nnzs.begin(); it != col_i_nnzs.end(); it++)
		{
			if (abs(col_i[*it]) > eps)
			{
				L.m_idx[k+1][i] = *it; //col k nonzero indices of L are stored
				L.m_x[k+1][i] = col_i[*it] / d; //col k nonzero values of L are stored
				L.list[*it].push_back(k+1); //update Llist
				count++;
				i++;
			}
		}
		// resize again, some 0s are dropped
		L.m_idx[k+1].resize(i);
		L.m_x[k+1].resize(i);

		count = count + 2;

		//advance list and L.first
		L.advance_first(k+1);
		advance_list(k+1);
		
		// ------------- reset vectors back to zero -----------------//
		for (auto it = col_i_nnzs.begin(); it != col_i_nnzs.end(); it++)
			col_i[*it] = 0;
		col_i_nnzs.clear();
		for (auto it = col_r_nnzs.begin(); it != col_r_nnzs.end(); it++)
			col_r[*it] = 0;
		col_r_nnzs.clear();
	}

	//assign number of non-zeros in L to L.nnz_count
	L.nnz_count = count;
}

#endif