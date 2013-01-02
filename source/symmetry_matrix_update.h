#ifndef SYMMETRY_MATRIX_UPDATE_H
#define SYMMETRY_MATRIX_UPDATE_H

using std::abs;


//----------------Column updates------------------//

template <class el_type>
inline void symmetry_matrix<el_type>::update_single(const int& j, const el_type& l_ki, const el_type& d, std::vector<el_type>& work, std::vector<int>& curr_nnzs, ultriangular_matrix<el_type>& L, vector<bool>& in_set) {
	//find where L(k, k+1:n) starts
	unsigned int i, offset = L.column_first[j];

	el_type factor = l_ki * d;
	for (i = offset; i < L.m_idx[j].size(); ++i) {
		work[L.m_idx[j][i]] -= factor * L.m_x[j][i];
	}
	
	//merge current non-zeros of col k with nonzeros of col *it.
	unordered_inplace_union(curr_nnzs, L.m_idx[j].begin() + offset,  L.m_idx[j].end(), in_set);
}

/*! \brief Performs a delayed update of subcolumn A(k:n,r). Result is stored in work vector. Nonzero elements of the work vector are stored in curr_nnzs.
	\param r the column number to be updated.
	\param work the vector for which all delayed-updates are computed to.
	\param curr_nnzs the nonzero elements of work.
	\param L the (partial) lower triangular factor of A.
	\param D the (partial) diagonal factor of A.
	\param in_set temporary storage for use in merging two lists of nonzero indices.
*/
template <class el_type>
inline void symmetry_matrix<el_type>::update(const int& r, std::vector<el_type>& work, std::vector<int>& curr_nnzs, ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, vector<bool>& in_set) {
	unsigned int j;
	int blk_sz;
	el_type d_12, l_ri;	

	//iterate across non-zeros of row k using Llist
	for (int i = 0; i < (int) L.list[r].size(); ++i) {
		j = L.list[r][i];
		
		l_ri = L.coeff(r, j, L.column_first[j]);
		update_single(j, l_ri, D[j], work, curr_nnzs, L, in_set); //update col using d11
		
		blk_sz = D.block_size(j);
		if (blk_sz == 2) {
			d_12 = D.off_diagonal(j);
			update_single(j + 1, l_ri, d_12, work, curr_nnzs, L, in_set);
		} else if (blk_sz == -2) {
			d_12 = D.off_diagonal(j-1);
			update_single(j - 1, l_ri, d_12, work, curr_nnzs, L, in_set); //update col using d12
		}
		
	}
}

#endif