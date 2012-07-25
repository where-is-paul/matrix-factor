template <class el_type>
inline void lilc_matrix<el_type> :: pivot(vector<idx_it>& swapk, vector<idx_it>& swapr, vector<list_it>& swapk_, vector<list_it>& swapr_, idx_vector_type& all_swaps, vector<bool>& in_set, elt_vector_type& col_k, idx_vector_type& col_k_nnzs, elt_vector_type& col_r, idx_vector_type& col_r_nnzs, lilc_matrix<el_type>& L, const int& k, const int& r)
{	
	//clear old stuff
	std::deque<int> row_k, row_r;
	std::pair<idx_it, elt_it> its_k, its_r;
	int i, j, idx, offset;
	
	row_k.clear();
	row_r.clear();
	
	col_k.clear();
	col_k_nnzs.clear();
	
	col_r.clear();
	col_r_nnzs.clear();
	//
	
	//----------pivot A ----------//
	swapr_.clear();
	swapk_.clear();
	all_swaps.clear();
	
	for (auto it = list[k].begin(); it != list[k].begin() + first[k]; it++) {
		row_r.push_back(*it);
	}
	
	for (auto it = list[r].begin(); it != list[r].begin() + first[r]; it++) {
		row_k.push_back(*it);
	}
	
	// cout << "first[k]: " << first[k] << " now outputting list[k]" << endl;
	// for (auto it = list[k].begin(); it != list[k].end(); it++) {
		// cout << *it << " ";
	// }
	// cout << endl;
	
	// cout << "first[r]: " << first[r] << " now outputting list[r]" << endl;
	// for (auto it = list[r].begin(); it != list[r].end(); it++) {
		// cout << *it << " ";
	// }
	// cout << endl;
	
	if (first[k] <= (int) list[k].size())
	all_swaps.assign(list[k].begin(), list[k].begin() + first[k]);

	if (first[r] <= (int) list[r].size())
	unordered_inplace_union(all_swaps, list[r].begin(), list[r].begin() + first[r], in_set);
	
	// cout << "sym swap on " << k << " and " << r << endl;
	for (auto it = all_swaps.begin(); it != all_swaps.end(); it++) {
		// cout << "curr swap on col " << *it << endl;
		
		// cout << m_idx[*it] << endl;
		safe_swap(m_idx[*it], k, r);
		// cout << m_idx[*it] << endl;
	}
	// cout << endl;
	all_swaps.clear();
	
	//after sym. perm, a_rr will be swapped to a_kk, so we put a_rr as first
	//elem of col k if its non-zero. this also means that we ensure the first
	//elem of col k is the diagonal element if it exists.
	if (coeff(r, r) !=0){
		col_k_nnzs.push_back(k);
		col_k.push_back(coeff(r, r));
	}
	
	//same as above, put a_kk in new col r if it exists.
	if (coeff(k, k) !=0){
		col_r_nnzs.push_back(r);
		col_r.push_back(coeff(k, k));
	}
	
	
	//first[r] should have A(r, k:r), not including A(k, r) element
	for (i = first[r]; i < (int) list[r].size(); i++) {
		if (coeffRef(r, list[r][i], its_k)) {
			if (list[r][i] == k) {
				col_k_nnzs.push_back(r);
			} else {
				col_k_nnzs.push_back(list[r][i]);
			}
			col_k.push_back(*its_k.second);
			
			*its_k.first = m_idx[list[r][i]].back();
			*its_k.second = m_x[list[r][i]].back();
			
			m_idx[list[r][i]].pop_back();
			m_x[list[r][i]].pop_back();
		}
	}

	if (m_idx[r].size() > 0) {
		offset = (m_idx[r][0] == r ? 1 : 0);
		std::copy(m_x[r].begin()+offset, m_x[r].end(), std::back_inserter(col_k));
		std::copy(m_idx[r].begin()+offset, m_idx[r].end(), std::back_inserter(col_k_nnzs));

		for (auto it = m_idx[r].begin() + offset; it != m_idx[r].end(); it++) {
			if (first[*it] < (int) list[*it].size())
			for (i = first[*it]; i < (int) list[*it].size(); i++) {
				if (list[*it][i] == r) {
					swapk_.push_back(list[*it].begin() + i);
					all_swaps.push_back(*it);
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
					ensure_invariant(idx, k, list[idx], first[idx], true);
					std::swap(list[idx][first[idx]], list[idx][list[idx].size() - 1]);
					list[idx].pop_back();
				}
				
				row_r.push_back(idx);
			} else if (idx > r) {
				col_r.push_back(m_x[k][i]);
				col_r_nnzs.push_back(idx);
				
				for (j = first[idx]; j < (int) list[idx].size(); j++) {
					if (list[idx][j] == k) {
						swapr_.push_back(list[idx].begin() + j);
						all_swaps.push_back(idx);
						break;
					}
				}
			}
		}
	}

	for (auto it = swapk_.begin(); it != swapk_.end(); it++) {
		**it = k;
	}
	
	for (auto it = swapr_.begin(); it != swapr_.end(); it++) {
		**it = r;
	}

	for (auto it = all_swaps.begin(); it != all_swaps.end(); it++) {
		ensure_invariant(*it, k, list[*it], first[*it], true);
	}

	for (auto it = col_k_nnzs.begin(); it != col_k_nnzs.end(); it++) {
		if ((*it != k) && (*it <= r)) {
			list[*it].push_back(k);
			ensure_invariant(*it, k, list[*it], first[*it], true);
		}
	}
	
	//set the kth col
	m_idx[k].swap(col_k_nnzs);
	m_x[k].swap(col_k);
	
	//set the rth col
	m_idx[r].swap(col_r_nnzs);
	m_x[r].swap(col_r);
	
	list[k].swap(row_k);
	list[r].swap(row_r);

	std::swap(first[k], first[r]);
	//--------end pivot A---------//
	
	//----------pivot L ----------//
	swapr.clear();
	swapk.clear();
	all_swaps.clear();
	
	for (auto it = L.list[k].begin(); it != L.list[k].end(); it++)
	{
		//L.ensure_invariant(*it, k, L.m_idx[*it], L.first[*it]);
		for (i = L.first[*it]; i < (int) L.m_idx[*it].size(); i++) {
			if (L.m_idx[*it][i] == k) {
				swapr.push_back(L.m_idx[*it].begin() + i);
				break;
			}
		}
	}
	
	for (auto it = L.list[r].begin(); it != L.list[r].end(); it++) {
		//L.ensure_invariant(*it, k, L.m_idx[*it], L.first[*it]);
		for (i = L.first[*it]; i < (int) L.m_idx[*it].size(); i++) {
			if (L.m_idx[*it][i] == r) {				
				swapk.push_back(L.m_idx[*it].begin() + i);
				break;
			}
		}
	}
	
	for (auto it = swapk.begin(); it != swapk.end(); it++) {
		**it = k;
	}
	
	for (auto it = swapr.begin(); it != swapr.end(); it++) {
		**it = r;
	}
	
	all_swaps.assign(L.list[r].begin(), L.list[r].end());
	unordered_inplace_union(all_swaps, L.list[k].begin(), L.list[k].end(), in_set);
	
	//invariant: ensure m_idx[:][L.first[i]+1] all contain the index nearest in value to k. 
	for (auto it = all_swaps.begin(); it != all_swaps.end(); it++) {
		L.ensure_invariant(*it, k, L.m_idx[*it], L.first[*it]);		
	}
	
	L.list[k].swap(L.list[r]);
	//--------end pivot L---------//
}
