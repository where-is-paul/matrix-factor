/*! \brief Computes the maximum (in absolute value) element of v(curr_nnzs) and it's index.
	\param v the vector whose max element is to be computed.
	\param curr_nnzs a list of indices representing non-zero elements in v.
	\param r the index of the maximum element of v
	\return the max element of v.
*/
template <class el_type>
inline double max(typename std::vector<el_type>& v, typename std::vector<int>& curr_nnzs, int& r) { 
	typename std::vector<int>::iterator it;
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
template <class el_type>
inline double norm(typename std::vector<el_type>& v, typename std::vector<int>& curr_nnzs) { 
	typename std::vector<int>::iterator it;
	double res = 0;
	for (it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) {
		res += pow(std::abs(v[*it]), 2);  
	}
	
	return sqrt(res);
}

/*! \brief Functor for comparing elements by value (in decreasing order) instead of by index.
	\param v the vector that contains the values being compared.
*/
template <class el_type>
struct by_value {
	std::vector<el_type>& v; 
	by_value(std::vector<el_type>& vec) : v(vec) {}
	bool operator()(int const &a, int const &b) const { 
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

template <class InputContainer, class InputIterator>
inline void unordered_inplace_union(InputContainer& a, InputIterator const& b_start, InputIterator const& b_end, vector<bool>& in_set)
{
	for (auto it = a.begin(); it != a.end(); it++) {
		in_set[*it] = true;
	}
	
	for (auto it = b_start; it != b_end; it++) {
		if (!in_set[*it]) {
			in_set[*it] = true;
			a.push_back(*it);
		}
	}
}

//-------------Dropping rules-------------//

/*! \brief Performs the dual-dropping criteria outlined in Li & Saad (2005).
	\param v the vector that whose elements will be selectively dropped.
	\param curr_nnzs the non-zeros in the vector v.
	\param lfil a parameter to control memory usage. Each column is guarannted to have fewer than lfil elements.
	\param tol a parameter to control agressiveness of dropping. Elements less than tol*norm(v) are dropped.
*/
template <class el_type>
inline void drop_tol(std::vector<el_type>& v, std::vector<int>& curr_nnzs, const int& lfil, const double& tol) { 
	typename std::vector<int>::const_iterator it;
	
	//determine dropping tolerance. all elements with value less than tolerance = tol * norm(v) is dropped.
	el_type tolerance = tol*norm<el_type>(v, curr_nnzs);
	for (it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) 
	if (std::abs(v[*it]) < tolerance) v[*it] = 0;
	
	//sort the remaining elements by value in decreasing order.
	by_value<el_type> sorter(v);
	std::sort(curr_nnzs.begin(), curr_nnzs.end(), sorter);
	
	//sort the first lfil elements by index, only these will be assigned into L.
	std::sort(curr_nnzs.begin(), curr_nnzs.begin() + std::min(lfil, (int) curr_nnzs.size()));
}

//----------------Column updates------------------//

template <class el_type>
inline void update_single(const int& k, const int& j, const el_type& l_ki, const el_type& d, std::vector<el_type>& work, std::vector<int>& curr_nnzs, lilc_matrix<el_type>& L, vector<bool>& in_set, bool include_kth = false) {
	//find where L(k, k+1:n) starts
	unsigned int i, offset = L.first[j];
	//if (offset >= L.m_idx[*it].size()) continue;
	if (L.m_idx[j][offset] < k) offset++;  //bug with L.first. shouldnt need more than one offset++.
	if (L.m_idx[j][offset] == k && !include_kth) offset++;
	
	for (i = offset; i < L.m_idx[j].size(); i++) {
		work[L.m_idx[j][i]] -= l_ki * d * L.m_x[j][i];
	}
	
	//merge current non-zeros of col k with nonzeros of col *it. 
	unordered_inplace_union(curr_nnzs, L.m_idx[j].begin() + offset,  L.m_idx[j].end(), in_set);
}

/*! \brief Performs a delayed update of subcolumn A(k+1:n,k). Result is stored in work vector. Nonzero elements of the work vector are stored in curr_nnzs.
	\param k the column number to be updated.
	\param work the vector for which all delayed-updates are computed to.
	\param curr_nnzs the nonzero elements of work.
	\param L the (partial) lower triangular factor of A.
	\param D the (partial) diagonal factor of A.
*/
template <class el_type>
inline void update(const int& k, std::vector<el_type>& work, std::vector<int>& curr_nnzs, lilc_matrix<el_type>& L, block_diag_matrix<el_type>& D, vector<bool>& in_set, bool include_kth = false) {
	unsigned int j;
	int blk_sz;
	el_type d_12, l_ki;	

	typename std::deque<int>::const_iterator it;
	//iterate across non-zeros of row k using Llist
	for (it = L.list[k].begin(); it != L.list[k].end(); it++) {
		j = *it;
		
		l_ki = L.coeff(k, j);
		update_single(k, j, l_ki, D[j], work, curr_nnzs, L, in_set, include_kth); //update col using d11
		
		blk_sz = D.block_size(j);
		if (blk_sz == 2) {
			d_12 = D.off_diagonal(j);
			update_single(k, j + 1, l_ki, d_12, work, curr_nnzs, L, in_set, include_kth);
		} else if (blk_sz == -2) {
			d_12 = D.off_diagonal(j-1);
			update_single(k, j - 1, l_ki, d_12, work, curr_nnzs, L, in_set, include_kth); //update col using d12
		}
		
	}
	
	for (auto it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) {
		in_set[*it] = false;
	}
}

//not needed anymore
template <class el_type>
inline void vec_add(std::vector<el_type>& v1, std::vector<int>& v1_nnzs, std::vector<el_type>& v2, std::vector<int>& v2_nnzs) {
	//merge current non-zeros of col k with nonzeros of col *it. 
	inplace_union(v1_nnzs, v2_nnzs.begin(), v2_nnzs.end());
	for (auto it = v1_nnzs.begin(); it != v1_nnzs.end(); it++) {
		v1[*it] += v2[*it];
	}
}

inline void safe_swap(std::vector<int>& curr_nnzs, int k, int r) {
	bool con_k = false, con_r = false;
	std::vector<int>::iterator k_idx, r_idx;
	for (auto it = curr_nnzs.begin(); it!= curr_nnzs.end(); it++) {
		if (*it == k) {
			con_k = true;
			k_idx = it;
		}
		
		if (*it == r) {
			con_r = true;
			r_idx = it;
		}
	}
	
	if (con_k == con_r) {
		//do nothing
	} else if (con_k) {
		*k_idx = r;
	} else {
		*r_idx = k;
	}
}
