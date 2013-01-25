#ifndef SYMMETRY_MATRIX_UPDATE_H
#define SYMMETRY_MATRIX_UPDATE_H

using std::abs;

//----------------Column updates------------------//

template <class el_type>
inline void symmetry_matrix<el_type>::update_single(const int& j, const el_type& l_rj, const el_type& d, std::vector<el_type>& work, std::vector<int>& work_nnzs, ultriangular_matrix<el_type>& L, vector<bool>& in_set)
{
	//find where L(k, k+1:n) starts
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

/*! \brief Performs a delayed update of subcolumn A(k:n,r). Result is stored in work vector. Nonzero elements of the work vector are stored in work_nnzs.
	\param r the column number to be updated.
	\param work the vector for which all delayed-updates are computed to.
	\param work_nnzs the nonzero elements of work.
	\param L the (partial) lower triangular factor of A.
	\param D the (partial) diagonal factor of A.
	\param in_set temporary storage for use in merging two lists of nonzero indices.
*/
template <class el_type>
inline void symmetry_matrix<el_type>::update(const int& r, std::vector<el_type>& work, std::vector<int>& work_nnzs, ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, vector<bool>& in_set)
{
	int j;
	int blk_sz;
	el_type d_12, l_rj;

	for (auto it = work_nnzs.begin(), end = work_nnzs.end(); it != end; it++)
		in_set[*it] = true;

	//iterate across non-zeros of row r using Llist
	for (int i = 0; i < (int) L.list[r].size(); i++) {
		j = L.list[r][i];
		
		l_rj = L.coeff(r, j, L.column_first[j]);
		update_single(j, l_rj, D[j], work, work_nnzs, L, in_set); //update col using d11
		
		blk_sz = D.block_size(j);
		if (blk_sz == 2) {
			d_12 = D.off_diagonal(j);
			update_single(j + 1, l_rj, d_12, work, work_nnzs, L, in_set);
		} else if (blk_sz == -2) {
			d_12 = D.off_diagonal(j-1);
			update_single(j - 1, l_rj, d_12, work, work_nnzs, L, in_set); //update col using d12
		}
	}

	for (auto it = work_nnzs.begin(), end = work_nnzs.end(); it != end; it++)
		in_set[*it] = false;
}

#endif