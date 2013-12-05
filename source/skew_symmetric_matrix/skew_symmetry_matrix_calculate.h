#ifndef SKEW_SYMMETRY_MATRIX_CALCULATE_H
#define SKEW_SYMMETRY_MATRIX_CALCULATE_H

using std::abs;

struct diagonal_filter {
  int k;
  diagonal_filter(int idx) : k(idx) {}
  bool operator()(const int& i) const {
    return (i == k) || (i == k+1);
  }
};

template <class el_type>
inline void skew_symmetry_matrix<el_type>::calculate(int& k, ultriangular_matrix<el_type>& L, skew_block_diag_matrix<el_type>& D, elt_vector_type& work1, elt_vector_type& work2, idx_vector_type& work1_nnzs, idx_vector_type& work2_nnzs, const double& tol, int& lfil, const double& stat_piv)
{
	//assign diagonal element to D
	D[k] = work1[k+1];
	work1[k] = work1[k+1] = work2[k] = work2[k+1] = 0;

	diagonal_filter filter(k);
	// erase diagonal, sub-diagonal, and super-diagonal elements from non-zero indices (to exclude it from being dropped)
	work1_nnzs.erase(std::remove_if(work1_nnzs.begin(), work1_nnzs.end(), filter), work1_nnzs.end());
	work2_nnzs.erase(std::remove_if(work2_nnzs.begin(), work2_nnzs.end(), filter), work2_nnzs.end());
	//performs the dual dropping procedure.
	drop_tol(work1, work1_nnzs, lfil, tol);
	drop_tol(work2, work2_nnzs, lfil, tol);
		
	if (abs(D[k]) < eps) D[k] = stat_piv; //statically pivot

	//resize kth column of L to proper size.
	L.m_x[k].resize(work2_nnzs.size());
	L.m_idx[k].assign(work2_nnzs.begin(), work2_nnzs.end());
	// fill L(k+2:n, k)
	int i = 0;
	for (idx_it it = work2_nnzs.begin(); it != work2_nnzs.end(); it++, i++)
	{
		L.m_x[k][i] = -work2[*it] / D[k]; //col k nonzero values of L are stored
		L.list[*it].push_back(k); //update Llist
	}
	L.nnz_count += work2_nnzs.size();
	

	//resize (k+1)th column of L to proper size.
	L.m_x[k+1].resize(work1_nnzs.size());
	L.m_idx[k+1].assign(work1_nnzs.begin(), work1_nnzs.end());	
	// fill L(k+2:n, k+1)
	i = 0;
	for (idx_it it = work1_nnzs.begin(); it != work1_nnzs.end(); it++, i++)
	{
		L.m_x[k+1][i] = work1[*it] / D[k]; //col k nonzero values of L are stored
		L.list[*it].push_back(k+1); //update Llist
	}
	L.nnz_count += work1_nnzs.size();
	
	//advance list and L.first
	L.advance_column(k+1);
	advance_list(k+1);
		
	// ------------- reset works back to zero -----------------//
	for (idx_it it = work1_nnzs.begin(); it != work1_nnzs.end(); it++)
		work1[*it] = 0;
	work1_nnzs.clear();
	for (idx_it it = work2_nnzs.begin(); it != work2_nnzs.end(); it++)
		work2[*it] = 0;
	work2_nnzs.clear();
}

#endif