using std::abs;

template <class el_type>
inline void lilc_matrix<el_type> :: pivot(swap_struct<el_type> s, vector<bool>& in_set, lilc_matrix<el_type>& L, const int& k, const int& r)
{	
	//clear old stuff
	std::vector<int> row_k, row_r;
	std::pair<idx_it, elt_it> its_k, its_r;
	int i, j, idx, offset;
	
	//for vectors of primitive types, clear is always constant time regardless of how many elements are in the container.
	s.col_clear();
	
	//----------pivot A ----------//
	s.swap_clear();
	
	for (auto it = list[k].begin(), end = list[k].begin() + first[k]; it != end; ++it) {
		row_r.push_back(*it);
	}
	
	for (auto it = list[r].begin(), end = list[r].begin() + first[r]; it != end; ++it) {
		row_k.push_back(*it);
	}
	
	s.all_swaps.assign(list[k].begin(), list[k].begin() + first[k]);

	unordered_inplace_union(s.all_swaps, list[r].begin(), list[r].begin() + first[r], in_set);
	
	for (auto it = s.all_swaps.begin(), end = s.all_swaps.end(); it != end; ++it) {
	
		safe_swap(m_idx[*it], k, r);

	}
	s.all_swaps.clear();
	
	//after sym. perm, a_rr will be swapped to a_kk, so we put a_rr as first
	//elem of col k if its non-zero. this also means that we ensure the first
	//elem of col k is the diagonal element if it exists.
	el_type elem = coeff(r, r);
	if (abs(elem) > eps){
		s.col_k_nnzs.push_back(k);
		s.col_k.push_back(elem);
	}
	
	//same as above, put a_kk in new col r if it exists.
	elem = coeff(k, k);
	if (abs(elem) > eps){
		s.col_r_nnzs.push_back(r);
		s.col_r.push_back(elem);
	}
	
	
	//first[r] should have A(r, k:r), not including A(k, r) element
	for (i = first[r]; i < (int) list[r].size(); i++) {
		j = list[r][i];
		if (coeffRef(r, j, its_k)) {
			if (j == k) {
				s.col_k_nnzs.push_back(r);
				row_r.push_back(k);
			} else {
				s.col_k_nnzs.push_back(j);
			}
			s.col_k.push_back(*its_k.second);
			
			*its_k.first = m_idx[j].back();
			*its_k.second = m_x[j].back();
			
			m_idx[j].pop_back();
			m_x[j].pop_back();
		}
	}

	if (m_idx[r].size() > 0) {
		offset = (m_idx[r][0] == r ? 1 : 0);
		std::copy(m_x[r].begin()+offset, m_x[r].end(), std::back_inserter(s.col_k));
		std::copy(m_idx[r].begin()+offset, m_idx[r].end(), std::back_inserter(s.col_k_nnzs));

		for (auto it = m_idx[r].begin() + offset; it != m_idx[r].end(); it++) {
			if (first[*it] < (int) list[*it].size())
			for (i = first[*it]; i < (int) list[*it].size(); i++) {
				if (list[*it][i] == r) {
					s.swapk.push_back(list[*it].begin() + i);
					s.all_swaps.push_back(*it);
					break;
				}
			}
		}
	}

	if (m_idx[k].size() > 0) {
		//swap A(k:r, k) with A(r, k:r);
		offset = (m_idx[k][0] == k ? 1 : 0);
		for (i = offset; i < (int) m_idx[k].size(); i++) {
			idx = m_idx[k][i];
			if (idx < r) {
				m_idx[idx].push_back(r);
				m_x[idx].push_back(m_x[k][i]);
				
				if (first[idx] < (int) list[idx].size()) {
					ensure_invariant(idx, k, list[idx], true);
					std::swap(list[idx][first[idx]], list[idx][list[idx].size() - 1]);
					list[idx].pop_back();
				}
				
				row_r.push_back(idx);
			} else if (idx > r) {
				s.col_r.push_back(m_x[k][i]);
				s.col_r_nnzs.push_back(idx);
				
				for (j = first[idx]; j < (int) list[idx].size(); j++) {
					if (list[idx][j] == k) {
						s.swapr.push_back(list[idx].begin() + j);
						s.all_swaps.push_back(idx);
						break;
					}
				}
			}
		}
	}

	for (auto it = s.swapk.begin(); it != s.swapk.end(); it++) {
		**it = k;
	}
	
	for (auto it = s.swapr.begin(); it != s.swapr.end(); it++) {
		**it = r;
	}

	for (auto it = s.all_swaps.begin(); it != s.all_swaps.end(); it++) {
		ensure_invariant(*it, k, list[*it], true);
	}

	for (auto it = s.col_k_nnzs.begin(); it != s.col_k_nnzs.end(); it++) {
		if ((*it != k) && (*it <= r)) {
			list[*it].push_back(k);
			ensure_invariant(*it, k, list[*it], true);
		}
	}
	
	//set the kth col
	m_idx[k].swap(s.col_k_nnzs);
	m_x[k].swap(s.col_k);
	
	//set the rth col
	m_idx[r].swap(s.col_r_nnzs);
	m_x[r].swap(s.col_r);
	
	list[k].swap(row_k);
	list[r].swap(row_r);

	std::swap(first[k], first[r]);
	//--------end pivot A---------//
	
	//----------pivot L ----------//
	s.swap_clear();

	for (auto it = L.list[k].begin(); it != L.list[k].end(); it++)
	{
		for (i = L.first[*it]; i < (int) L.m_idx[*it].size(); i++) {
			if (L.m_idx[*it][i] == k) {
				s.swapr.push_back(L.m_idx[*it].begin() + i);
				break;
			}
		}
	}
	
	for (auto it = L.list[r].begin(); it != L.list[r].end(); it++) {
		for (i = L.first[*it]; i < (int) L.m_idx[*it].size(); i++) {
			if (L.m_idx[*it][i] == r) {				
				s.swapk.push_back(L.m_idx[*it].begin() + i);
				break;
			}
		}
	}
	
	for (auto it = s.swapk.begin(); it != s.swapk.end(); it++) {
		**it = k;
	}
	
	for (auto it = s.swapr.begin(); it != s.swapr.end(); it++) {
		**it = r;
	}
	
	s.all_swaps.assign(L.list[r].begin(), L.list[r].end());
	unordered_inplace_union(s.all_swaps, L.list[k].begin(), L.list[k].end(), in_set);
	
	for (auto it = s.all_swaps.begin(); it != s.all_swaps.end(); it++) {
		L.ensure_invariant(*it, k, L.m_idx[*it]);		
	}
	
	L.list[k].swap(L.list[r]);
	//--------end pivot L---------//
}
