#ifndef _CSC_MATRIX_ILDL_H_
#define _CSC_MATRIX_ILDL_H_

#include <algorithm>
#include <cmath>
#include <list>
#include <deque>

/*! possible optimizations/to-do's: 
	- some easy parallelization of for loops. 
	
	- using a vector to accumulate all nonzeros in previous cols
	  during iteration.
	  
	- small changes to loops (e.g. assign .begin(), and .end() so that they will not be evaluated each iter.
	---> for (auto itr = new_values.begin(), end_itr = new_values.end(); itr != end_itr; ++itr ) )
	
	- make a special diagonal matrix class for storing 1x1 and 2x2 pivots later.
	
*/

/*! \brief Computes the norm of v(curr_nnzs)
	\param v the vector whose norm is to be computed.
	\param curr_nnzs a list of indices representing non-zero elements in v.
	\return the norm of v.
*/
template <class idx_type, class el_type>
inline double norm(typename std::vector<el_type>& v, typename std::vector<idx_type>& curr_nnzs) { 
    typename std::vector<idx_type>::iterator it;
    double res = 0;
    for (it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) {
        res += pow(std::abs(v[*it]), 2);  
    }
    
    return res;
}

/*! \brief Functor for comparing elements by value (in decreasing order) instead of by index.
	\param v the vector that contains the values being compared.
*/
template <class idx_type, class el_type>
struct by_value {
    std::vector<el_type>& v; 
    by_value(std::vector<el_type>& vec) : v(vec) {}
    bool operator()(idx_type const &a, idx_type const &b) const { 
        return std::abs(v[a]) > std::abs(v[b]);
    }
};

/*! \brief Performs the dual-dropping criteria outlined in Li & Saad (2005).
	\param v the vector that whose elements will be selectively dropped.
	\param curr_nnzs the non-zeros in the vector v.
	\param lfil a parameter to control memory usage. Each column is guarannted to have fewer than lfil elements.
	\param tol a parameter to control agressiveness of dropping. Elements less than tol*norm(v) are dropped.
*/
template <class idx_type, class el_type>
inline void drop_tol(std::vector<el_type>& v, std::vector<idx_type>& curr_nnzs, const int& lfil, const double& tol) { 
    typename std::vector<idx_type>::iterator it;
	
	//determine dropping tolerance. all elements with value less than tolerance = tol * norm(v) is dropped.
    el_type tolerance = tol*norm<idx_type, el_type>(v, curr_nnzs);
    for (it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) 
        if (std::abs(v[*it]) < tolerance) v[*it] = 0;
     
	//sort the remaining elements by value in decreasing order.
    by_value<idx_type, el_type> sorter(v);
    std::sort(curr_nnzs.begin(), curr_nnzs.end(), sorter);
    
	//sort the first lfil elements by index, only these will be assigned into L.
    std::sort(curr_nnzs.begin(), curr_nnzs.begin() + std::min(lfil, (int) curr_nnzs.size()));
}

/*! \brief Performs an inplace union of two sorted lists (a and b), removing duplicates in the final list.
	\param a the sorted list to contain the final merged list.
	\param b_start an iterator to the start of b.
	\param b_end an iterator to the end of b.
*/
template <class InputContainer, class InputIterator>
inline void inplace_union(InputContainer& a, InputIterator const& b_start, InputIterator const& b_end)
{
    int mid = a.size(); //store the end of first sorted range

    //copy the second sorted range into the destination vector
    std::copy(b_start, b_end, std::back_inserter(a));

    //perform the in place merge on the two sub-sorted ranges.
    std::inplace_merge(a.begin(), a.begin() + mid, a.end());

    //remove duplicate elements from the sorted vector
    a.erase(std::unique(a.begin(), a.end()), a.end());
}


template <class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: ildl(csc_matrix<idx_type, el_type>& L, elt_vector_type& D, int lfil, double tol)
{	
	//Llist is a deque of linked lists that gives the non-zero elements in each row of L. since at any time we may swap between two rows, we require a linked lists for each row of L. a deque is used as it might be desireable to deallocate all linked lists for rows i < k on step k (this is currently not done, as the memory used in maintaining linked lists for all rows is not much).
	std::vector< std::deque< idx_type > > Llist; 
	
	//work is a work vector for the current column. Lfirst is a linked list that gives the first nonzero element in column k with row index i > k. (i.e. the first nonzero in L(k+1:n, k).
	std::vector<el_type> work(n_cols(), 0), Lfirst(n_cols(), 0);
	
	std::vector<idx_type> curr_nnzs; //non-zeros on current col.
	
	int count = 0; //the current non-zero in L.
	idx_type i, j, k, offset;
	el_type l_ik;
	typename std::deque<idx_type>::const_iterator it;
	
	curr_nnzs.reserve(n_cols()); //makes sure that there is enough space if every element in the column is nonzero
	Llist.resize(n_cols()); //allocate a vector of size n for Llist.
	
	L.resize(n_rows(), n_cols(), (lfil+1)*n_cols()); //(+1 because there are 1s on the diagonal. they wont need to be stored if we want to optimize)
	D.resize(n_cols()); 
	
	for (i = 0; i < n_cols(); i++) {
		D[i] = coeff(i,i);
	}
	
	for (k = 0; k < n_cols(); k++) {
	    //zero out work vector
	    std::fill (work.begin() + k, work.end(), 0);
	    
	    //the +1 avoids assigning diagonal element as nonzero since its stored in D. the min is in case m_col_idx[k] == m_col_idx[k+1] (an empty column).
	    curr_nnzs.assign (m_row_idx.begin() + min(m_col_idx[k] + 1, m_col_idx[k+1]), m_row_idx.begin() + m_col_idx[k+1]);
		
		//assigns the non zeros in A(k,:) to the work vector. since only the lower diagonal of A is stored, this is essentially A(k,k+1:n).
	    for (j = m_col_idx[k] + 1; j < m_col_idx[k+1]; j++) {
			work[m_row_idx[j]] = m_x[j];
		}
		
		//iterate across non-zeros of row k using Llist
		if (k < n_cols() - 1)
		for (it = Llist[k].begin(); it != Llist[k].end(); it++) { //Llist[k] contains the row idx of row k.
			
			//find where L(k, k+1:n) starts
			offset = L.m_col_idx[*it] + Lfirst[*it] + 1; 
			if (L.m_row_idx[offset] == k) {
				Lfirst[*it]++; //update Lfirst
			
				for (j = offset+1; j < L.m_col_idx[*it+1]; j++) {
					if (L.m_row_idx[j] > k)
						work[L.m_row_idx[j]] -= L.m_x[offset] * D[*it] * L.m_x[j];
				}
				
				
				//find exactly where to merge with binary search.
				typename std::vector<idx_type>::iterator start = upper_bound(L.m_row_idx.begin() + L.m_col_idx[*it] + 1, L.m_row_idx.begin() + L.m_col_idx[*it+1], k);
				
				//merge current non-zeros of col k with nonzeros of col *it. 
				inplace_union(curr_nnzs, start, L.m_row_idx.begin() + L.m_col_idx[*it+1]);
				
			}
			
		}		
		
		//performs the dual dropping procedure.
		drop_tol(work, curr_nnzs, lfil, tol);	
		
		//get 1s on the diagonal
		L.m_row_idx[count] = k;
		L.m_x[count] = 1;
		count++;
		
		
		if (k < n_cols() - 1)
		for (i = 0; i < (idx_type) std::min(lfil, (int) curr_nnzs.size()); i++) {
		    L.m_row_idx[count] = curr_nnzs[i]; //row_idx of L is updated
		    L.m_x[count] = work[curr_nnzs[i]]/D[k]; //work vector is scaled by D[k]
			Llist[curr_nnzs[i]].push_back(k); //update Llist
			count++;
		}
		
		L.m_col_idx[k+1] = count; //the end of the current column is assigned to col_idx
		
		if (k < n_cols() - 1) {
			//finds out where L(k+1:n, k) starts
			offset = L.m_col_idx[k] + Lfirst[k] + 1;
			for (i = offset; i < L.m_col_idx[k+1]; i++) {
				l_ik = L.m_x[i];
				D[L.m_row_idx[i]] -= l_ik * D[k] * l_ik;	//update diagonal
			}
		}
	}
	
	//resize vectors of L down to the right size.
	L.m_row_idx.resize(count); 
	L.m_x.resize(count);
	
}

#endif 
