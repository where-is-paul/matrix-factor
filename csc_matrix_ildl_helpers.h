/*! \brief Computes the maximum (in absolute value) element of v(curr_nnzs) and it's index.
	\param v the vector whose max element is to be computed.
	\param curr_nnzs a list of indices representing non-zero elements in v.
	\param r the index of the maximum element of v
	\return the max element of v.
*/
template <class idx_type, class el_type>
inline double max(typename std::vector<el_type>& v, typename std::vector<idx_type>& curr_nnzs, idx_type& r) { 
    typename std::vector<idx_type>::iterator it;
    double res = 0;
    for (it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) {
		if (std::abs(v[*it]) > res) {
			res = std::abs(v[*it]);
			r = *it;
		}
    }
	
    return res;
}

/*! \brief Computes the norm of v(curr_nnzs).
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

//-------------Dropping rules-------------//

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

//----------------Column updates------------------//


/*! \brief Performs a delayed update of subcolumn A(k+1:n,k). Result is stored in work vector. Nonzero elements of the work vector are stored in curr_nnzs.
	\param k the column number to be updated.
	\param work the vector for which all delayed-updates are computed to.
	\param curr_nnzs the nonzero elements of work.
	\param L the (partial) lower triangular factor of A.
	\param D the (partial) diagonal factor of A.
	\param Lfirst vector containing location of first nonzero of column k that is below the diagonal.
	\param Llist vector containing nonzeros in each row of L.
*/
template <class idx_type, class el_type>
inline void update(const idx_type& k, std::vector<el_type>& work, std::vector<idx_type>& curr_nnzs, csc_matrix<idx_type, el_type>& L, std::vector<el_type>& D, std::vector<idx_type>& Lfirst, std::vector< std::deque< idx_type > >& Llist) {
	idx_type j, offset;
	typename std::deque<idx_type>::const_iterator it;
	//iterate across non-zeros of row k using Llist
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
}
