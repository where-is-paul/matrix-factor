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
inline void drop_tol(std::vector<el_type>& v, const idx_type& k, std::vector<idx_type>& curr_nnzs, const int& lfil, const double& tol) { 
    //optimize later
    typename std::vector<el_type>::iterator it;
    el_type tolerance = tol*norm<idx_type, el_type>(v.begin()+k, v.end());
    for (it = v.begin()+k; it != v.end(); it++) 
        if (std::abs(*it) < tolerance) *it = 0;
        
    //do lfil later. approach: use merge while doing delayed updates to obtain list of non zeros. sort curr_nnzs by value. only assign lfil largest to new vector
    by_value<idx_type, el_type> sorter(v);
    std::sort(curr_nnzs.begin(), curr_nnzs.end(), sorter);
    
    std::sort(curr_nnzs.begin(), curr_nnzs.begin() + std::min(lfil, (int) curr_nnzs.size()));
}

template <class InputContainer, class InputIterator>
inline void inplace_union(InputContainer& a, InputIterator b_start, InputIterator b_end)
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
	std::vector<idx_type> curr_nnzs;
	std::deque<idx_type> all_nnzs;
	
	int count = 0; 
	idx_type i, j, k, offset;
	el_type l_ik;
	typename std::deque<idx_type>::const_iterator it;
	
	curr_nnzs.reserve(n_cols());
	
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
	    curr_nnzs.assign (m_row_idx.begin() + min(m_col_idx[k] + 1, m_col_idx[k+1]), m_row_idx.begin() + m_col_idx[k+1]); //the min is in case col_idx[k] == col_idx[k+1]
		
	    //needs optimization for sparse vector add later (maybe turn work vector into sparse vector,
	    //use STL merge on row indices).
		//#pragma omp parallel for
	    for (j = m_col_idx[k] + 1; j < m_col_idx[k+1]; j++) {
			work[m_row_idx[j]] = m_x[j];
		}
		
		//get 1s on the diagonal
		L.m_row_idx[count] = k;
		L.m_x[count] = 1;
		count++;
		
		if (k < n_cols() - 1) {
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
					
				}
				
			}	
			
			//merge non-zero pattern of all accumulated previous cols with curr col. binary search to find starting point
			inplace_union(curr_nnzs, upper_bound(all_nnzs.begin(), all_nnzs.end(), k), all_nnzs.end() );
			
			//<--- drop the elements that are <= k in index! (implement linked list later), when optimizing, no need to erase and reallocate every iter.
			curr_nnzs.erase( remove_if(curr_nnzs.begin(), curr_nnzs.end(), [&k] (idx_type i) -> bool {return i <= k;}), curr_nnzs.end() );	
			
			drop_tol(work, k, curr_nnzs, lfil, tol);	

			for (i = 0; i < (idx_type) std::min(lfil, (int) curr_nnzs.size()); i++) {
				L.m_row_idx[count] = curr_nnzs[i];
				L.m_x[count] = work[curr_nnzs[i]]/D[k];
				
				Llist[curr_nnzs[i]].push_back(k); //update Llist
				
				count++;
			}
		
			//offset does not need to be calculated. its really just count (before count was changed by the for loop above.)
			offset = L.m_col_idx[k] + Lfirst[k] + 1;
			
			for (i = offset; i < L.m_col_idx[k+1]; i++) {	//use Lfirst here later.
				l_ik = L.m_x[i];
				D[L.m_row_idx[i]] -= l_ik * D[k] * l_ik;	
			}
		}
		
		L.m_col_idx[k+1] = count;
		
		//make more explicit to the fact that it is accumulating unions. perhaps shorten all_nnzs as we go by removing elements <= k. have to change to deque then.
		inplace_union(all_nnzs, L.m_row_idx.begin() + L.m_col_idx[k] + 1, L.m_row_idx.begin() + L.m_col_idx[k+1]);
		
		if (all_nnzs[0] == k) all_nnzs.pop_front();
	}
	
	L.m_row_idx.resize(count);
	L.m_x.resize(count);
	
}

#endif 
