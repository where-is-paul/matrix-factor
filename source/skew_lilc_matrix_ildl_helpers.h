#ifndef _SKEW_LILC_MATRIX_ILDL_HELPERS_H
#define _SKEW_LILC_MATRIX_ILDL_HELPERS_H

using std::abs;

//----------------Column updates------------------//

template <class el_type>
inline void skew_update_single(const int& r, const int& i, const el_type& l_ri, const el_type& d, std::vector<el_type>& work, std::vector<int>& work_nnzs, lilc_matrix<el_type>& L, vector<bool>& in_set)
{
	unsigned int j;
	el_type factor = l_ri * d;
	for (j = L.first[i]; j < L.m_idx[i].size(); j++)
		work[L.m_idx[i][j]] -= factor * L.m_x[i][j];
	
	//merge current non-zeros of col k with nonzeros of col *it.
	unordered_inplace_union(work_nnzs, L.m_idx[i].begin() + L.first[i],  L.m_idx[i].end(), in_set);
}

/*! \brief Performs a delayed update of subcolumn A(k:n,r). Result is stored in work vector. Nonzero elements of the work vector are stored in curr_nnzs.
	\A(k:n, r) += L(k:n, i+1)*D(i, i+1)*L(r, i) + L(k:n, i)*D(i+1, i)*L(r, i+1) (i = 0, 2, ..., k-2)
	\param r the column number to be updated.
	\param work the vector for which all delayed-updates are computed to.
	\param work_nnzs the nonzero elements of work.
	\param L the (partial) lower triangular factor of A.
	\param D the (partial) diagonal factor of A.
	\param in_set temporary storage for use in merging two lists of nonzero indices.
*/
template <class el_type>
inline void skew_update(const int& r, std::vector<el_type>& work, std::vector<int>& work_nnzs, lilc_matrix<el_type>& L, skew_block_diag_matrix<el_type>& D, vector<bool>& in_set)
{
	unsigned int i;
	el_type l_ri;

	//iterate across non-zeros of row k using Llist
	for (int j = 0; j < (int) L.list[r].size(); j++)
	{
		i = L.list[r][j];
		
		l_ri = L.coeff(r, i, L.first[i]);
		if (i%2 == 0)
			skew_update_single(r, i+1, l_ri, D[i], work, work_nnzs, L, in_set);
		else
			skew_update_single(r, i-1, l_ri, -D[i-1], work, work_nnzs, L, in_set);
	}
}

#endif