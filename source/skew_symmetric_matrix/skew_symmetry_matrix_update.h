#ifndef SKEW_SYMMETRY_MATRIX_UPDATE_H
#define SKEW_SYMMETRY_MATRIX_UPDATE_H

using std::abs;

//----------------Column updates------------------//

template <class el_type>
inline void skew_symmetry_matrix<el_type>::update_single(const int& j, const el_type& l_rj, const el_type& d, std::vector<el_type>& work, std::vector<int>& work_nnzs, ultriangular_matrix<el_type>& L, vector<bool>& in_set)
{

	el_type factor = l_rj * d;
	for (int i = L.column_first[j]; i < (int) L.m_idx[j].size(); i++)
	{
		work[L.m_idx[j][i]] -= factor * L.m_x[j][i];
		//merge current non-zeros of col k with nonzeros of col *it.
		if (!in_set[L.m_idx[j][i]])
		{
			in_set[L.m_idx[j][i]] = true;
			work_nnzs.push_back(L.m_idx[j][i]);
		}
	}
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
inline void skew_symmetry_matrix<el_type>::update(const int& r, std::vector<el_type>& work, std::vector<int>& work_nnzs, ultriangular_matrix<el_type>& L, skew_block_diag_matrix<el_type>& D, vector<bool>& in_set)
{
	int j;
	el_type l_rj;

	for (idx_it it = work_nnzs.begin(), end = work_nnzs.end(); it != end; it++)
		in_set[*it] = true;

	//iterate across non-zeros of row r using Llist
	for (int i = 0; i < (int) L.list[r].size(); i++)
	{
		j = L.list[r][i];
		
		l_rj = L.coeff(r, j, L.column_first[j]);
		if (j%2 == 0)
			update_single(j+1, l_rj, D[j], work, work_nnzs, L, in_set);
		else
			update_single(j-1, l_rj, -D[j-1], work, work_nnzs, L, in_set);
	}

	for (idx_it it = work_nnzs.begin(), end = work_nnzs.end(); it != end; it++)
		in_set[*it] = false;
}

#endif