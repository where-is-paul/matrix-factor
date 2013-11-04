#ifndef SYMMETRY_MATRIX_CALCULATE_H
#define SYMMETRY_MATRIX_CALCULATE_H

using std::abs;

struct diagonal_filter {
  int k;
  diagonal_filter(int idx) : k(idx) {}
  bool operator()(const int& i) const {
    return (i == k) || (i == k+1);
  }
};

template <class el_type>
inline void symmetry_matrix<el_type>::calculate(int& k, ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, elt_vector_type& col_i, elt_vector_type& col_r, idx_vector_type& col_i_nnzs, idx_vector_type& col_r_nnzs, bool& size_two_piv, const double& tol, int& lfil, vector<bool>& in_set, const double& stat_piv)
{
	el_type det_D, D_inv11, D_inv22, D_inv12;	//for use in 2x2 pivots
	el_type l_11, l_12;							//for use in 2x2 pivots
	
	//assign diagonal element to D
	D[k] = col_i[k];
	//performs the dual dropping procedure.
	if (!size_two_piv)
	{
		//erase diagonal element from non-zero indices (to exclude it from being dropped)
		col_i_nnzs.erase(std::remove(col_i_nnzs.begin(), col_i_nnzs.end(), k), col_i_nnzs.end());
		//perform dual dropping criteria on col_i
		drop_tol(col_i, col_i_nnzs, lfil, tol);
	}
	else
	{
    diagonal_filter filter(k);
		//erase diagonal 2x2 block from non-zero indices (to exclude it from being dropped)
		col_i_nnzs.erase(std::remove_if(col_i_nnzs.begin(), col_i_nnzs.end(), filter), col_i_nnzs.end());
		col_r_nnzs.erase(std::remove_if(col_r_nnzs.begin(), col_r_nnzs.end(), filter), col_r_nnzs.end());
		//assign pivot to D
		D.off_diagonal(k) = col_i[k+1];
		D[k+1] = col_r[k+1];

		//compute inverse of the 2x2 block diagonal pivot.
		det_D = col_i[k] * col_r[k+1] - col_i[k+1] * col_i[k+1];
		if (abs(det_D) < eps) {
			det_D = stat_piv;  //statically pivot;
		}
		D_inv11 = col_r[k+1] / det_D;
		D_inv22 = col_i[k] / det_D;
		D_inv12 = -col_i[k+1] / det_D;
		
		//merge nonzeros of col_i_nnzs and col_r_nnzs together so iterating through them will be easier
		unordered_inplace_union(col_i_nnzs, col_r_nnzs.begin(), col_r_nnzs.end(), in_set);

		//multiply inverse of pivot to col_i and col_r (gives us two columns of l)
		for (idx_it it = col_i_nnzs.begin(); it != col_i_nnzs.end(); it++)
		{
			l_11 = col_i[*it]*D_inv11 + col_r[*it]*D_inv12;
			l_12 = col_i[*it]*D_inv12 + col_r[*it]*D_inv22;
				
			//note that col_i and col_r roughly share the same non-zero indices
			col_i[*it] = l_11;
			col_r[*it] = l_12;
		}
		
		//since the col_i and col_r non-zero indices are roughly the same,
		//we can copy it over to col_r_nnzs
		col_r_nnzs.assign(col_i_nnzs.begin(), col_i_nnzs.end());
		
		//perform dual dropping procedure on col_i and col_r
		drop_tol(col_i, col_i_nnzs, lfil, tol);
		drop_tol(col_r, col_r_nnzs, lfil, tol);
	}

	//resize kth column of L to proper size.
	L.m_idx[k].resize(col_i_nnzs.size());
	L.m_x[k].resize(col_i_nnzs.size());

	if (!size_two_piv)
	{
		if (abs(D[k]) < eps) {
			D[k] = stat_piv; //statically pivot
		}
		int i = 0;
		for (idx_it it = col_i_nnzs.begin(); it != col_i_nnzs.end(); it++, i++)
		{ 
			L.m_idx[k][i] = *it; //col k nonzero indices of L are stored
			L.m_x[k][i] = col_i[*it] / D[k]; //col k nonzero values of L are stored
			L.list[*it].push_back(k); //update Llist
		}
		L.nnz_count += col_i_nnzs.size();

		//advance list and L.first
		L.advance_column(k);
		advance_list(k);
	}
	else
	{
		//resize k+1th column of L to proper size.
		L.m_idx[k+1].resize(col_r_nnzs.size());
		L.m_x[k+1].resize(col_r_nnzs.size());

		int i = 0;
		for (idx_it it = col_i_nnzs.begin(); it != col_i_nnzs.end(); it++, i++)
		{
			L.m_idx[k][i] = *it; //col k nonzero values of L are stored
			L.m_x[k][i] = col_i[*it]; //col k nonzero indices of L are stored
			L.list[*it].push_back(k); //update L.list
		}
			
		i = 0;
		for (idx_it it = col_r_nnzs.begin(); it != col_r_nnzs.end(); it++, i++)
		{
			L.m_idx[k+1][i] = *it; //col k+1 nonzero values of L are stored
			L.m_x[k+1][i] = col_r[*it]; //col k+1 nonzero indices of L are stored
			L.list[*it].push_back(k+1); //update L.list
		}
		L.nnz_count += col_i_nnzs.size() + col_r_nnzs.size();

		//update list and L.first
		L.advance_column(k+1);
		advance_list(k+1);
	}
	
	// ------------- reset col_i and col_r back to zero -----------------//
	for (idx_it it = col_i_nnzs.begin(); it != col_i_nnzs.end(); it++)
		col_i[*it] = 0;
	col_i_nnzs.clear(); //zero out col_i vector
	for (idx_it it = col_r_nnzs.begin(); it != col_r_nnzs.end(); it++)
		col_r[*it] = 0;
	col_r_nnzs.clear(); //zero out col_r vector
	
	if (size_two_piv)
	{
		k++;
		size_two_piv = false;
	}
}

#endif