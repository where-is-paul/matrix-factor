#ifndef _CSC_MATRIX_ILDL_H_
#define _CSC_MATRIX_ILDL_H_

#include <algorithm>
#include <cmath>
#include <list>
#include <deque>

template <class idx_type, class el_type>
inline double norm(typename std::vector<el_type>::iterator v_start, typename std::vector<el_type>::iterator v_end) { 
    //optimize later
    typename std::vector<el_type>::iterator it;
    double res = 0;
    //may want to use abs later if want to deal with complex
	//#pragma omp parallel for reduction(+:res)
    for (it = v_start; it != v_end; it++) {
        res += pow(std::abs(*it), 2);  
    }
    
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
inline void drop_tol(std::vector<el_type>& v, const idx_type& k, std::vector<idx_type>& nnzs, const int& lfil, const double& tol) { 
    //optimize later
    typename std::vector<el_type>::iterator it;
    el_type tolerance = tol*norm<idx_type, el_type>(v.begin()+k, v.end());
    for (it = v.begin()+k; it != v.end(); it++) 
        if (std::abs(*it) < tolerance) *it = 0;
        
    //do lfil later. approach: use merge while doing delayed updates to obtain list of non zeros. sort nnzs by value. only assign lfil largest to new vector
    by_value<idx_type, el_type> sorter(v);
    std::sort(nnzs.begin(), nnzs.end(), sorter);
    
    std::sort(nnzs.begin(), nnzs.begin() + std::min(lfil, (int) nnzs.size()));
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
	std::vector< std::deque< idx_type > > Llist;
	std::vector<el_type> work(n_cols(), 0), Lfirst(n_cols(), 0);
	std::vector<idx_type> nnzs;
	
	int count = 0; 
	idx_type i, j, k, offset;
	el_type l_ik;
	typename std::deque<idx_type>::const_iterator it;
	
	nnzs.reserve(n_cols());
	
	Llist.resize(n_cols()); //reserve nspots for Llist.
	
	L.resize(n_rows(), n_cols(), (lfil+1)*n_cols()); //(+1 because there are 1s on the diagonal. they wont need to be stored if we want to optimize)
	D.resize(n_cols()); //make a tridiagonal matrix class later? will need for 2x2 pivots.
	
	//#pragma omp parallel for
	for (i = 0; i < n_cols(); i++) {
		D[i] = coeff(i,i);
	}
	
	for (k = 0; k < n_cols(); k++) {
	    //zero out work vector (opt. later)
	    std::fill (work.begin() + k, work.end(), 0); //can change to work.begin() + k, but need to fix the norm after.
	    
		
	    //the +1 is to avoid assigning diagonal element as nonzero since its stored in D already
	    nnzs.assign (m_row_idx.begin() + min(m_col_idx[k] + 1, m_col_idx[k+1]), m_row_idx.begin() + m_col_idx[k+1]); //the min is in case col_idx[k] == col_idx[k+1]
		
	    //needs optimization for sparse vector add later (maybe turn work vector into sparse vector,
	    //use STL merge on row indices).
		//#pragma omp parallel for
	    for (j = m_col_idx[k] + 1; j < m_col_idx[k+1]; j++) {
			work[m_row_idx[j]] = m_x[j];
		}
		
		if (k < n_cols() - 1)
		for (it = Llist[k].begin(); it != Llist[k].end(); it++) {
			//update Lfirst
			offset = L.m_col_idx[*it] + Lfirst[*it] + 1;
			
			if (L.m_row_idx[offset] == k) {
				Lfirst[*it]++; //update Lfirst
			
				//assumes matrix is symmetric with only lower triangular part stored. may want to change later.
				for (j = offset+1; j < L.m_col_idx[*it+1]; j++) {
					if (L.m_row_idx[j] > k)
						work[L.m_row_idx[j]] -= L.m_x[offset] * D[*it] * L.m_x[j];
				}
				
				
				//find exactly where to merge with binary search. this function will probably be customized to deal with the linked lists
				inplace_union(nnzs, L.m_row_idx.begin() + L.m_col_idx[*it] + 1, L.m_row_idx.begin() + L.m_col_idx[*it+1]);
				
			}
			
		}	
		
		//<--- drop the elements that are <= k in index! (implement linked list later), when optimizing, no need to erase and reallocate every iter.
		nnzs.erase( remove_if(nnzs.begin(), nnzs.end(), [&k] (idx_type i) -> bool {return i <= k;}), nnzs.end() );	
		
		drop_tol(work, k, nnzs, lfil, tol);	
		
		//get 1s on the diagonal
		L.m_row_idx[count] = k;
		L.m_x[count] = 1;
		count++;
		
		
		if (k < n_cols() - 1)
		for (i = 0; i < (idx_type) std::min(lfil, (int) nnzs.size()); i++) {
		    L.m_row_idx[count] = nnzs[i];
		    L.m_x[count] = work[nnzs[i]]/D[k];
			
			Llist[nnzs[i]].push_back(k); //update Llist
			
			count++;
		}
		
		L.m_col_idx[k+1] = count;
		
		if (k < n_cols() - 1) {
			offset = L.m_col_idx[k] + Lfirst[k] + 1;
			
			for (i = offset; i < L.m_col_idx[k+1]; i++) {	//use Lfirst here later.
				l_ik = L.m_x[i];
				D[L.m_row_idx[i]] -= l_ik * D[k] * l_ik;	
			}
		}
		
		//put Lfirst and Llist updates here
	}
	
	L.m_row_idx.resize(count);
	L.m_x.resize(count);
	
}

#endif 
