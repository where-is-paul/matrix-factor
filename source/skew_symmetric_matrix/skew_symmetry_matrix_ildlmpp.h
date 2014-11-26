#ifndef SKEW_SYMMETRY_MATRIX_ILDLMPP_H
#define SKEW_SYMMETRY_MATRIX_ILDLMPP_H

using std::abs;
using std::cout;
using std::min;
using std::endl;

template <class el_type>
void skew_symmetry_matrix<el_type>::ildlmpp(ultriangular_matrix<el_type>& L, skew_block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol, const double& pp_tol)
{
	//----------------- initialize temporary variables --------------------//
	const int ncols = n_cols(); // number of cols in A.
	const double stat_piv = 1e-6;

	int lfil;
	if (fill_factor >= 1e4) lfil = ncols; // just incase users decide to enter a giant fill factor for fun...
	else lfil = (int) (fill_factor*nnz()/ncols); // roughly a factor of 2 since only lower tri. of A is stored

	el_type w1, w2;
	bool firstcol = true;
	int j, k, r, r1, r2;

	vector<bool> in_set(ncols, false); // bitset used for unsorted merges
	square_matrix_swap_struct<el_type> s;	// struct containing temp vars used in pivoting.

	elt_vector_type work1(ncols, 0), work2(ncols, 0); // work vector for the current columns
	idx_vector_type work1_nnzs, work2_nnzs;  // non-zeros on current col.
	work1_nnzs.reserve(ncols); // reserves space for worse case (entire col is non-zero)
	work2_nnzs.reserve(ncols);

	//--------------- allocate memory for L and D ------------------//
	L.resize(ncols, ncols, min(2*lfil, 32)); //allocate a vector of size n for Llist as well
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
		update(k, work1, work1_nnzs, L, D, in_set);
		update(k+1, work2, work2_nnzs, L, D, in_set);
		
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
			L.advance_column(k);
			// swap rows and columns of k+1 and r
			if (k+1 != r)
			{
				for (idx_it it = work2_nnzs.begin(); it != work2_nnzs.end(); it++)
					work2[*it] = 0;
				work2_nnzs.clear();
				for (j = list_first[r]; j < (int) list[r].size(); j++)
				{
					work2_nnzs.push_back(list[r][j]);
					work2[list[r][j]] = -this->coeff(r, list[r][j]);
				}
				work2_nnzs.insert(work2_nnzs.end(), m_idx[r].begin(), m_idx[r].end());
				for (j = 0; j < (int) m_idx[r].size(); j++)
					work2[m_idx[r][j]] = m_x[r][j];

				update(r, work2, work2_nnzs, L, D, in_set);
				pivot(s, L, k+1, r);
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
			for (idx_it it = work1_nnzs.begin(); it != work1_nnzs.end(); it++)
				work1[*it] = 0;
			work1_nnzs.clear();
			for (j = list_first[r]; j < (int) list[r].size(); j++)
			{
				work1_nnzs.push_back(list[r][j]);
				work1[list[r][j]] = -this->coeff(r, list[r][j]);
			}
			work1_nnzs.insert(work1_nnzs.end(), m_idx[r].begin(), m_idx[r].end());
			for (j = 0; j < (int) m_idx[r].size(); j++)
				work1[m_idx[r][j]] = m_x[r][j];

			update(r, work1, work1_nnzs, L, D, in_set);
			pivot(s, L, k, r);
			std::swap(perm[k], perm[r]);
			std::swap(work1[k], work1[r]);
			std::swap(work2[k], work2[r]);
			safe_swap(work1_nnzs, k, r);
			safe_swap(work2_nnzs, k, r);
			advance_list(k);
			L.advance_column(k);
		}
		//--------------end pivoting--------------//

		calculate(k, L, D, work1, work2, work1_nnzs, work2_nnzs, tol, lfil, stat_piv);
	}

	L.nnz_count += ncols; // diagonal 1s are not stored in L, but L.nnz_count should include them
}

#endif