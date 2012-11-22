#ifndef _SKEW_LILC_MATRIX_ILDLMPP_H_
#define _SKEW_LILC_MATRIX_ILDLMPP_H_

template <class el_type>
void skew_lilc_matrix<el_type>::ildlmpp(lilc_matrix<el_type>& L, skew_block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol)
{
	//----------------- initialize temporary variables --------------------//
	const int ncols = n_cols(); // number of cols in A.
	int lfil = (int)(2 * fill_factor*nnz() / ncols); // roughly a factor of 2 since only lower tri. of A is stored
	el_type w1, w2, d;
	bool firstcol = true;
	int count = 0; // the total number of nonzeros stored in L.
	int i, j, k, r, r1, r2;
	vector<bool> in_set(ncols, false); // bitset used for unsorted merges
	skew_swap_struct<el_type> s;	// struct containing temp vars used in pivoting.
	elt_vector_type work1(ncols, 0), work2(ncols, 0); // work vector for the current columns
	idx_vector_type work1_nnzs, work2_nnzs;  // non-zeros on current col.
	
	work1_nnzs.reserve(ncols); // reserves space for worse case (entire col is non-zero)
	work2_nnzs.reserve(ncols);

	//--------------- allocate memory for L and D ------------------//
	L.resize(ncols, ncols);
	D.resize(ncols );

	//------------------- main loop: factoring begins -------------------------//
	for (k = 0; k < ncols; k=k+2)
	{
		// assign nonzeros indices and values of A(k+1:n, k) to work1_nnzs
		work1_nnzs.assign(m_idx[k].begin(), m_idx[k].end());
		for (j = 0; j < (int) work1_nnzs.size(); j++)
			work1[work1_nnzs[j]] = m_x[k][j];

		// assign nonzeros indices and values of A(k+2:n, k+1) to work2_nnzs
		work2_nnzs.assign(m_idx[k+1].begin(), m_idx[k+1].end());
		for (j = 0; j < (int) work2_nnzs.size(); j++)
			work2[work2_nnzs[j]] = m_x[k+1][j];
		if (abs(work1[k+1]) > eps) // add the superdiagonal element to work2, no need to add the value
			work2_nnzs.push_back(k);

		// do delayed updates on current column
		// A(k:n, r) += L(k:n, i+1)*D(i, i+1)*L(r, i) + L(k:n, i)*D(i+1, i)*L(r, i+1) (i = 0, 2, ..., k-2)
		skew_update(k, work1, work1_nnzs, L, D, in_set);
		skew_update(k+1, work2, work2_nnzs, L, D, in_set);
		
		work2[k] = -work1[k+1]; // the correct superdiagonal value is set here
		
		//find maximum element in work and store its index in r.
		w1 = max(work1, work1_nnzs, r1);
		w2 = max(work2, work2_nnzs, r2);
		
		if ((w1 >= w2) || (r2 == k))
		{
			r = r1;
			firstcol = true;
		}
		else
		{
			r = r2;
			firstcol = false;
		}

		//--------------begin pivoting--------------//
		if (w1 < eps && w2 < eps)
			cout << "singular" << endl;
		else if (firstcol) // maximum element in col k
		{
			advance_list(k);
			L.advance_first(k);
			// swap rows and columns of k+1 and r
			if (k+1 != r)
			{
				for (auto it = work2_nnzs.begin(); it != work2_nnzs.end(); it++)
					work2[*it] = 0;
				work2_nnzs.clear();
				for (j = first[r]; j < (int) list[r].size(); j++)
				{
					work2_nnzs.push_back(list[r][j]);
					work2[list[r][j]] = -coeff(r, list[r][j]);
				}
				work2_nnzs.insert(work2_nnzs.end(), m_idx[r].begin(), m_idx[r].end());
				for (j = 0; j < (int) m_idx[r].size(); j++)
					work2[m_idx[r][j]] = m_x[r][j];

				skew_update(r, work2, work2_nnzs, L, D, in_set);
				pivot(s, in_set, L, k+1, r);
				std::swap(perm[k+1], perm[r]);
				std::swap(work1[k+1], work1[r]);
				std::swap(work2[k+1], work2[r]);
				safe_swap(work1_nnzs, k+1, r);
				safe_swap(work2_nnzs, k+1, r);
			}
		}
		else // maximum element in col k+1
		{
			// swap rows and columns of k and r
			for (auto it = work1_nnzs.begin(); it != work1_nnzs.end(); it++)
				work1[*it] = 0;
			work1_nnzs.clear();
			for (j = first[r]; j < (int) list[r].size(); j++)
			{
				work1_nnzs.push_back(list[r][j]);
				work1[list[r][j]] = -coeff(r, list[r][j]);
			}
			work1_nnzs.insert(work1_nnzs.end(), m_idx[r].begin(), m_idx[r].end());
			for (j = 0; j < (int) m_idx[r].size(); j++)
				work1[m_idx[r][j]] = m_x[r][j];

			skew_update(r, work1, work1_nnzs, L, D, in_set);
			pivot(s, in_set, L, k, r);
			std::swap(perm[k], perm[r]);
			std::swap(work1[k], work1[r]);
			std::swap(work2[k], work2[r]);
			safe_swap(work1_nnzs, k, r);
			safe_swap(work2_nnzs, k, r);
			advance_list(k);
			L.advance_first(k);
		}
		//--------------end pivoting--------------//
		D[k] = d = work1[k+1];
		work1[k] = work1[k+1] = work2[k] = work2[k+1] = 0;

		// erase diagonal, sub-diagonal, and super-diagonal elements from non-zero indices (to exclude it from being dropped)
		work1_nnzs.erase(std::remove_if(work1_nnzs.begin(), work1_nnzs.end(), isKorKp1(k)), work1_nnzs.end());
		work2_nnzs.erase(std::remove_if(work2_nnzs.begin(), work2_nnzs.end(), isKorKp1(k)), work2_nnzs.end());
		//performs the dual dropping procedure.
		drop_tol(work1, work1_nnzs, lfil, tol);
		drop_tol(work2, work2_nnzs, lfil, tol);
		
		//resize kth column of L to proper size.
		L.m_idx[k].resize(work2_nnzs.size()+1);
		L.m_x[k].resize(work2_nnzs.size()+1);
		L.m_idx[k][0] = k;
		L.m_x[k][0] = 1;
		
		if (abs(D[k]) < eps) D[k] = 1e-6; //statically pivot
		// fill L(k+2:n, k)
		i = 1;
		for (auto it = work2_nnzs.begin(); it != work2_nnzs.end(); it++)
		{
			if (abs(work2[*it]) > eps)
			{
				L.m_idx[k][i] = *it; //col k nonzero indices of L are stored
				L.m_x[k][i] = -work2[*it] / d; //col k nonzero values of L are stored
				L.list[*it].push_back(k); //update Llist
				count++;
				i++;
			}
		}
		// resize again, some 0s are dropped
		L.m_idx[k].resize(i);
		L.m_x[k].resize(i);


		//resize (k+1)th column of L to proper size.
		L.m_idx[k+1].resize(work1_nnzs.size()+1);
		L.m_x[k+1].resize(work1_nnzs.size()+1);
		L.m_idx[k+1][0] = k+1;
		L.m_x[k+1][0] = 1;
		
		// fill L(k+2:n, k+1)
		i = 1;
		for (auto it = work1_nnzs.begin(); it != work1_nnzs.end(); it++)
		{
			if (abs(work1[*it]) > eps)
			{
				L.m_idx[k+1][i] = *it; //col k nonzero indices of L are stored
				L.m_x[k+1][i] = work1[*it] / d; //col k nonzero values of L are stored
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
		
		// ------------- reset works back to zero -----------------//
		for (auto it = work1_nnzs.begin(); it != work1_nnzs.end(); it++)
			work1[*it] = 0;
		work1_nnzs.clear();
		for (auto it = work2_nnzs.begin(); it != work2_nnzs.end(); it++)
			work2[*it] = 0;
		work2_nnzs.clear();
	}

	//assign number of non-zeros in L to L.nnz_count
	L.nnz_count = count;
}

#endif