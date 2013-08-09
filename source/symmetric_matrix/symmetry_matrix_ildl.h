#ifndef SYMMETRY_MATRIX_ILDL_H
#define SYMMETRY_MATRIX_ILDL_H

using std::abs;

template <class el_type>
void symmetry_matrix<el_type>::ildl(ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol, const double& pp_tol)
{
	//----------------- initialize temporary variables --------------------//
	const int ncols = n_cols(); //number of cols in A.
	const double stat_piv = 1e-6;

	int lfil;
	if (fill_factor >= 1e4) lfil = ncols; //just incase users decide to enter a giant fill factor for fun...
	else lfil = (int) (fill_factor*nnz()/ncols); //roughly a factor of 2 since only lower tri. of A is stored
	
	const el_type alpha = (1.0+sqrt(17.0))/8.0;  //for use in pivoting.
	el_type wi(-1), wr(-1), d(-1);

	vector<bool> in_set(ncols, false); //bitset used for unsorted merges
	square_matrix_swap_struct<el_type> s;	//struct containing temp vars used in pivoting.
	
	elt_vector_type col_i(ncols, 0), col_r(ncols, 0); ////work vector for the current column
	idx_vector_type col_i_nnzs, col_r_nnzs;  //non-zeros on current col.
	col_i_nnzs.reserve(ncols); //reserves space for worse case (entire col is non-zero)
	col_r_nnzs.reserve(ncols); //reserves space for worse case (entire col is non-zero)

	int i, k, r;
	bool size_two_piv = false;	//boolean indicating if the pivot is 2x2 or 1x1
	
	//--------------- allocate memory for L and D ------------------//
	L.resize(ncols, ncols, 2*lfil); //allocate a vector of size n for Llist as well
	D.resize(ncols );
	

	//------------------- main loop: factoring begins -------------------------//
	for (k = 0; k < ncols; k++)
	{
		//assign nonzeros indices of A(k:n, k) to col_i_nnzs
		col_i_nnzs.assign(m_idx[k].begin(), m_idx[k].end());
		//assign nonzero values of A(k:n, k) to col_i
		for (i = 0; i < (int) col_i_nnzs.size(); i++)
			col_i[col_i_nnzs[i]] = m_x[k][i];

		//--------------begin pivoting--------------//

		//do delayed updates on current column. col_i = Sum_{i=0}^{k-1} L(k,i) * D(i,i) * L(k:n, i)
		//(the formula above generalizes to block matrix form in the case of 2x2 pivots).
		update(k, col_i, col_i_nnzs, L, D, in_set);
		
		//store diagonal element in di. set diagonal element in col_i vector to 0
		//since we want to find the maximum off-diagonal element.
		d = col_i[k]; col_i[k] = 0;
		//find maximum element in col_i and store its index in r.
		wi = max(col_i, col_i_nnzs, r);
    
    //we do partial pivoting here, where we take the first element u in the column that satisfies
    //|u| > pp_tol*|wi|. for more information, consult "A Partial Pivoting Strategy for Sparse 
    //Symmetric Matrix Decomposition" by J.H. Liu (1987).
    int t = r; //stores location of u 
    el_type u = wi; //stores value of u
    for (i = 0; i < (int) col_i_nnzs.size(); i++) {
      if (abs(col_i[col_i_nnzs[i]])-pp_tol*wi > eps ) {
        t = col_i_nnzs[i];
        u = col_i[t];
        break;
      }
    }
		col_i[k] = d;

		//bunch-kaufman partial pivoting is used below. for a more detailed reference,
		//refer to "Accuracy and Stability of Numerical Algorithms." by Higham (2002).
		//------------------- begin bunch-kaufman pivoting ------------------//
		if (wi < eps) {
			//case 0: do nothing. pivot is k.
		} else if ( (alpha * wi - abs(d)) < eps  ) {
			//case 1: do nothing. pivot is k.
		}
		else
		{
      //since we are doing partial pivoting, we should treat u and t like wi and r, so
      //we'll just reassign wi and r. note: this has to go in the else clause since
      //we still use the old wi for case 0 and case 1.
      wi = u;
      r = t;
      
			//assign all nonzero indices and values in A(r, k:r)
			//( not including A(r,r) ) to col_r and col_r_nnzs
			for (i = list_first[r]; i < (int) list[r].size(); i++) {
				col_r_nnzs.push_back(list[r][i]);
				col_r[list[r][i]] = this->coeff(r, list[r][i]);
			}

			//assign nonzero indices of A(r:n, r) to col_r_nnzs 
			col_r_nnzs.insert(col_r_nnzs.end(), m_idx[r].begin(), m_idx[r].end());
			//assign nonzero values of to col_r
			for (i = 0; i < (int) m_idx[r].size(); i++)
				col_r[m_idx[r][i]] = m_x[r][i];

			//perform delayed updates on col_r. col_r = Sum_{i=0}^{k-1} L(r,i) * D(i,i) * L(k:n, i).
			//(the formula above generalizes to block matrix form in the case of 2x2 pivots).
			update(r, col_r, col_r_nnzs, L, D, in_set);

			d = col_r[r]; col_r[r] = 0;
			//find maximum element in col_r.
			wr = max(col_r, col_r_nnzs, i);
			col_r[r] = d;

			if ((alpha*wi*wi - abs(col_i[k])*wr) < eps) {
				//case 2: do nothing. pivot is k.
			}
			else if ((alpha * wr - abs(d)) < eps)
			{
				//case 3: pivot is k with r: 1x1 pivot case.
				//--------pivot A and L ---------//
				this->pivot(s, L, k, r);

				//----------pivot rest ----------//

				//permute perm
				std::swap(perm[k], perm[r]);

				col_i.swap(col_r);	//swap col_i with col_r.
				std::swap(col_i[k], col_i[r]); //swap kth and rth row of col_i
				col_i_nnzs.swap(col_r_nnzs);	//swap col_i_nnzs with col_r_nnzs
				safe_swap(col_i_nnzs, k, r); //swap k and r if they are present in col_i_nnzs
				//--------end pivot rest---------//
			}
			else
			{
				//case 4: pivot is k+1 with r: 2x2 pivot case.
				//must advance list for 2x2 pivot since we are pivoting on col k+1
				this->advance_list(k);
				//for the same reason as above, we must advance L.first as well
				L.advance_column(k);

				size_two_piv = true;

				if (k+1 != r)
				{
					//symmetrically permute row/col k+1 and r.
					this->pivot(s, L, k+1, r);

					//----------pivot rest ----------//

					//permute perm
					std::swap(perm[k+1], perm[r]);

					//swap rows k+1 and r of col_i and col_r
					std::swap(col_i[k+1], col_i[r]);
					std::swap(col_r[k+1], col_r[r]);
					//swap k+1 and r in col_i_nnzs and col_r_nnzs
					safe_swap(col_i_nnzs, k+1, r);
					safe_swap(col_r_nnzs, k+1, r);
				}
			}
		}
		//--------------end pivoting--------------//

		calculate(k, L, D, col_i, col_r, col_i_nnzs, col_r_nnzs, size_two_piv, tol, lfil, in_set, stat_piv);
	}

	L.nnz_count += ncols; // diagonal 1s are not stored in L, but L.nnz_count should include them
}

#endif