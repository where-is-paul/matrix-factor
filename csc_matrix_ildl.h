#ifndef _CSC_MATRIX_ILDL_H_
#define _CSC_MATRIX_ILDL_H_

#include <algorithm>
#include <cmath>
#include <list>

template <class idx_type, class el_type>
inline double norm(std::vector<el_type>& v) { 
    //optimize later
    double res = 0;
    //may want to use abs later if want to deal with complex
	//#pragma omp parallel for reduction(+:res)
    for (idx_type i = 0; i < (idx_type) v.size(); i++) res += pow(std::abs(v[i]), 2);  
    return res;
}

//optimize later.
template <class idx_type, class el_type>
struct by_value {
    std::vector<el_type>& v; 
    by_value(std::vector<el_type>& vec) : v(vec) {}
    bool operator()(idx_type const &a, idx_type const &b) const { 
        return std::abs(v[a]) > std::abs(v[b]);
    }
};

template <class idx_type, class el_type>
inline void drop_tol(std::vector<el_type>& v, std::vector<idx_type>& nnzs, const int& lfil, const double& tol) { 
    //optimize later
    el_type tolerance = tol*norm<idx_type, el_type>(v);
    for (idx_type i = 0; i < (idx_type) v.size(); i++) 
        if (std::abs(v[i]) < tolerance) v[i] = 0;
        
    //do lfil later. approach: use merge while doing delayed updates to obtain list of non zeros. sort nnzs by value. only assign lfil largest to new vector
    by_value<idx_type, el_type> sorter(v);
    std::sort(nnzs.begin(), nnzs.end(), sorter);
    
    nnzs.resize( std::min(lfil, (int) nnzs.size()) );
    std::sort(nnzs.begin(), nnzs.end());
}

template <class idx_type, class InputIterator>
inline void inplace_union(std::vector<idx_type>& a, InputIterator b_start, InputIterator b_end)
{
    int mid = a.size(); //Store the end of first sorted range

    //First copy the second sorted range into the destination vector
    std::copy(b_start, b_end, std::back_inserter(a));

    //Then perform the in place merge on the two sub-sorted ranges.
    std::inplace_merge(a.begin(), a.begin() + mid, a.end());

    //Remove duplicate elements from the sorted vector
    a.erase(std::unique(a.begin(), a.end()), a.end());
}


template <class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: ildl(csc_matrix<idx_type, el_type>& L, elt_vector_type& D, int lfil, double tol)
{
	std::vector< std::list< idx_type > > Llist;
	std::vector<el_type> work(n_cols(), 0), Lfirst(n_cols(), 0);
	std::vector<idx_type> nnzs;
	
	int count = 0; idx_type i, j, k;
	nnzs.reserve(n_cols());
	
	L.resize(n_rows(), n_cols(), lfil*n_cols());
	D.resize(n_cols()); //make a tridiagonal matrix class later? will need for 2x2 pivots.
	
	//set d_k = a_kk switch to binary search later
	//#pragma omp parallel for
	for (i = 0; i < n_cols(); i++) {
		D[i] = coeff(i,i);
	}
	
	for (k = 0; k < n_cols(); k++) {
	    //zero out work vector (opt. later)
	    std::fill (work.begin() + k, work.end(), 0); //can change to work.begin() + k, but need to fix the norm after.
	    
	    //the +1 is to avoid assigning diagonal element as nonzero since its stored in D already
	    nnzs.assign (m_row_idx.begin() + m_col_idx[k] + 1, m_row_idx.begin() + m_col_idx[k+1]);
	    
		
	    //needs optimization for sparse vector add later (maybe turn work vector into sparse vector,
	    //use STL merge on row indices).
		//#pragma omp parallel for
	    for (j = m_col_idx[k]; j < m_col_idx[k+1]; j++) {
			if (m_row_idx[j] > k)
				work[m_row_idx[j]] = m_x[j];
		}
		
		
		if (k < n_cols() - 1)
		for (i = 0; i < k; i++) {
			//update Lfirst
			idx_type offset = L.m_col_idx[i] + Lfirst[i] + 1;
			if (L.m_row_idx[offset] == k) {
				Lfirst[i]++; //update Lfirst
			
				//assumes matrix is symmetric with only lower triangular part stored. may want to change later.
				for (j = offset+1; j < L.m_col_idx[i+1]; j++) {
					if (L.m_row_idx[j] > k)
						work[L.m_row_idx[j]] -= L.m_x[offset] * D[i] * L.m_x[j];
				}
				
				//find exactly where to merge with binary search. this function will probably be customized to deal with the linked lists
				inplace_union(nnzs, L.m_row_idx.begin() + L.m_col_idx[i] + 1, L.m_row_idx.begin() + L.m_col_idx[i+1]);
			}
			
		}	
		
		//<--- drop the elements that are <= k in index! (implement linked list later), when optimizing, no need to erase and reallocate every iter.
		nnzs.erase( remove_if(nnzs.begin(), nnzs.end(), [&k] (idx_type i) -> bool {return i <= k;}), nnzs.end() );	
		
		drop_tol(work, nnzs, lfil, tol);	
		
		//get 1s on the diagonal
		L.m_row_idx[count] = k;
		L.m_x[count] = 1;
		count++;
		
		if (k < n_cols() - 1)
		for (i = 0; i < (idx_type) nnzs.size(); i++) {
		    L.m_row_idx[count] = nnzs[i];
		    L.m_x[count] = work[nnzs[i]]/D[k];
		    count++;
		}
		
		L.m_col_idx[k+1] = count;
		
		
		el_type l_ik;
		for (i = k+1; i < n_cols(); i++) {	//use Lfirst here later.
			l_ik = L.coeff(i, k);
			if (l_ik != 0)
				D[i] -= l_ik * D[k] * l_ik;			
		}
		
		//put Lfirst and Llist updates here
	}
	
	L.m_row_idx.resize(count);
	L.m_x.resize(count);
	
}

#endif 
