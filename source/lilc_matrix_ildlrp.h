#ifndef _LILC_MATRIX_ILDLRP_H_
#define _LILC_MATRIX_ILDLRP_H_

template <class el_type>
void lilc_matrix<el_type>::ildlrp(lilc_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol)
{
	//----------------- initialize temporary variables --------------------//
	const int ncols = n_cols(); // number of cols in A.
	const double stat_piv = 1e-6;

	int lfil;
	if (fill_factor > 1e4) lfil = ncols; //just incase users decide to enter a giant fill factor for fun...
	else lfil = (int) (2*fill_factor*nnz()/ncols); //roughly a factor of 2 since only lower tri. of A is stored
	
	const el_type alpha = (1.0+sqrt(17.0))/8.0;  //for use in pivoting.
	el_type wi(-1), wr(-1), di(-1), dr(-1);
	el_type det_D, D_inv11, D_inv22, D_inv12;	//for use in 2x2 pivots
	el_type l_11, l_12;							//for use in 2x2 pivots

	vector<bool> in_set(ncols, false); //bitset used for unsorted merges
	swap_struct<el_type> s;	//struct containing col_r vars used in pivoting.
	elt_vector_type col_i(ncols, 0), col_r(ncols, 0); ////work vector for the current column
	idx_vector_type col_i_nnzs, col_r_nnzs;  //non-zeros on current col.
	col_i_nnzs.reserve(ncols); //reserves space for worse case (entire col is non-zero)
	col_r_nnzs.reserve(ncols);

	int count = 0; //the total number of nonzeros stored in L.
	int i, j, k, r, newr;
	bool size_two_piv = false;	//boolean indicating if the pivot is 2x2 or 1x1

	//--------------- allocate memory for L and D ------------------//
	L.resize(ncols, ncols);
	D.resize(ncols );
	
	//------------------- main loop: factoring begins -------------------------//
	for (k = 0; k < ncols; k++)
	{
		i = k;
		//assign nonzeros indices of A(k:n, i) to col_i_nnzs
		col_i_nnzs.assign(m_idx[k].begin(), m_idx[k].end());
		for (j = 0; j < (int) col_i_nnzs.size(); j++)
			col_i[col_i_nnzs[j]] = m_x[k][j];
		//do delayed updates on current column. work = Sum_{i=0}^{k-1} L(k,i) * D(i,i) * L(k:n, i)
		//(the formula above generalizes to block matrix form in the case of 2x2 pivots).
		update(i, col_i, col_i_nnzs, L, D, in_set);
		di = col_i[i]; col_i[i] = 0;
		wi = max(col_i, col_i_nnzs, r);

		if (alpha * wi <= abs(di) + eps)
			size_two_piv = false;
		else
		{
			while (true)
			{
				// assign nonzeros indices and values of A(k:n, r) to col_r_nnzs
				for (auto it = col_r_nnzs.begin(); it != col_r_nnzs.end(); it++)
					col_r[*it] = 0;
				col_r_nnzs.clear();
				for (j = first[r]; j < (int) list[r].size(); j++)
				{
					col_r_nnzs.push_back(list[r][j]);
					col_r[list[r][j]] = coeff(r, list[r][j]);
				}
				col_r_nnzs.insert(col_r_nnzs.end(), m_idx[r].begin(), m_idx[r].end());
				for (j = 0; j < (int) m_idx[r].size(); j++)
					col_r[m_idx[r][j]] = m_x[r][j];
				update(r, col_r, col_r_nnzs, L, D, in_set);

				//find maximum element in col_r and store its index in newr.
				dr = col_r[r]; col_r[r] = 0;
				wr = max(col_r, col_r_nnzs, newr);

				if (alpha * wr <= abs(dr) + eps)
				{
					size_two_piv = false;
					// swap rows and columns k and r
					pivot(s, in_set, L, k, r);
					std::swap(perm[k], perm[r]);
					std::swap(col_r[k], col_r[r]);
					safe_swap(col_r_nnzs, k, r);
					col_i.swap(col_r);
					col_i_nnzs.swap(col_r_nnzs);
					di = dr;
					break;
				}
				else if (abs(wi - wr) < eps)
				{
					size_two_piv = true;
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
					di = dr;
					wi = wr;
					r = newr;
					col_i.swap(col_r);
					col_i_nnzs.swap(col_r_nnzs);
				}
			}
		}
		//--------------end pivoting--------------//

		//erase diagonal element from non-zero indices (to exclude it from being dropped)
		col_i_nnzs.erase(std::remove(col_i_nnzs.begin(), col_i_nnzs.end(), k), col_i_nnzs.end());
		D[k] = di;

		//performs the dual dropping procedure.
		if (!size_two_piv)
			drop_tol(col_i, col_i_nnzs, lfil, tol);
		else
		{
			//erase diagonal 2x2 block from non-zero indices (to exclude it from being dropped)
			col_r_nnzs.erase(std::remove(col_r_nnzs.begin(), col_r_nnzs.end(), k), col_r_nnzs.end());
			col_i_nnzs.erase(std::remove(col_i_nnzs.begin(), col_i_nnzs.end(), k+1), col_i_nnzs.end());
			col_r_nnzs.erase(std::remove(col_r_nnzs.begin(), col_r_nnzs.end(), k+1), col_r_nnzs.end());
			//assign pivot to D
			D.off_diagonal(k) = col_i[k+1]; col_i[k+1] = 0; col_r[k] = 0;
			D[k+1] = dr;

			//compute inverse of the 2x2 block diagonal pivot.
			det_D = di * dr - D.off_diagonal(k) * D.off_diagonal(k);
			if (abs(det_D) < eps) det_D = stat_piv;  //statically pivot;
			D_inv11 = dr / det_D;
			D_inv22 = di / det_D;
			D_inv12 = -D.off_diagonal(k) / det_D;
			
			//merge nonzeros of curr and temp together so iterating through them will be easier
			unordered_inplace_union(col_i_nnzs, col_r_nnzs.begin(), col_r_nnzs.end(), in_set);

			//multiply inverse of pivot to work and temp (gives us two columns of l)
			for (auto it = col_i_nnzs.begin(); it != col_i_nnzs.end(); it++)
			{
				l_11 = col_i[*it]*D_inv11 + col_r[*it]*D_inv12;
				l_12 = col_i[*it]*D_inv12 + col_r[*it]*D_inv22;
				
				//note that col_i and col_r roughly share the same non-zero indices
				col_i[*it] = l_11;
				col_r[*it] = l_12;
			}
			
			//since the work and temp non-zero indices are roughly the same,
			//we can copy it over to temp_nnzs
			col_r_nnzs.assign(col_i_nnzs.begin(), col_i_nnzs.end());
			
			// perform dual dropping procedure
			drop_tol(col_i, col_i_nnzs, lfil, tol);
			drop_tol(col_r, col_r_nnzs, lfil, tol);
		}

		//resize kth column of L to proper size.
		L.m_idx[k].resize(col_i_nnzs.size()+1);
		L.m_x[k].resize(col_i_nnzs.size()+1);
		
		// assign 1s to diagonal of L.
		L.m_idx[k][0] = k;
		L.m_x[k][0] = 1;
		count++;
		
		if (!size_two_piv)
		{
			if ( abs(D[k]) < eps) D[k] = stat_piv; //statically pivot
			i = 1;
			for (auto it = col_i_nnzs.begin(); it != col_i_nnzs.end(); it++)
			{ 
				if (abs(col_i[*it]) > eps)
				{
					L.m_idx[k][i] = *it; //col k nonzero indices of L are stored
					L.m_x[k][i] = col_i[*it] / D[k]; //col k nonzero values of L are stored
					L.list[*it].push_back(k); //update Llist
					count++;
					i++;
				}
			}
			
			//advance list and L.first
			L.advance_first(k);
			advance_list(k);
		}
		else
		{
			// resize k+1th column of L to proper size.
			L.m_idx[k+1].resize(col_r_nnzs.size() + 1);
			L.m_x[k+1].resize(col_r_nnzs.size() + 1);

			//assign 1s to diagonal of L.
			L.m_idx[k+1][0] = k+1;
			L.m_x[k+1][0] = 1;
			count++;

			i = 1;
			for (auto it = col_i_nnzs.begin(); it != col_i_nnzs.end(); it++)
			{
				if (abs(col_i[*it]) > eps)
				{
					L.m_idx[k][i] = *it; //col k nonzero values of L are stored
					L.m_x[k][i] = col_i[*it]; //col k nonzero indices of L are stored
					L.list[*it].push_back(k); //update L.list
					count++;
					i++;
				}
			}
			
			j = 1;
			for (auto it = col_r_nnzs.begin(); it != col_r_nnzs.end(); it++)
			{
				if (abs(col_r[*it]) > eps)
				{
					L.m_idx[k+1][j] = *it; //col k+1 nonzero values of L are stored
					L.m_x[k+1][j] = col_r[*it]; //col k+1 nonzero indices of L are stored
					L.list[*it].push_back(k+1); //update L.list
					count++;
					j++;
				}
			}

			//update list and L.first
			L.advance_first(k+1);
			advance_list(k+1);
		}
		
		//resize columns of L to correct size
		L.m_idx[k].resize(i);
		L.m_x[k].resize(i);

		if (size_two_piv)
		{
			L.m_idx[k+1].resize(j);
			L.m_x[k+1].resize(j);
			k++;
		}

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