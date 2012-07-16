template <class el_type>
inline void lilc_matrix<el_type> :: pivot(vector<idx_it>& swapk, vector<idx_it>& swapr, vector<list_it>& swapk_, vector<list_it>& swapr_, idx_vector_type& all_swaps, vector<bool>& in_set, elt_vector_type& col_k, idx_vector_type& col_k_nnzs, elt_vector_type& col_r, idx_vector_type& col_r_nnzs, lilc_matrix<el_type>& L, const int& k, const int& r)
{	
	std::pair<idx_it, elt_it> its_k, its_r;
	int i, j, offset;
	
	col_k.clear();
	col_k_nnzs.clear();
	
	col_r.clear();
	col_r_nnzs.clear();
	
	//----------pivot A ----------//
	swapr_.clear();
	swapk_.clear();
	all_swaps.clear();
	
	//swap A(:, k) with A(:, r)
	//invariant: ensure m_idx[:][0] and m_x[:][0] contains the smallest index element.
	if (coeff(r, r) !=0){
		col_k_nnzs.push_back(k);
		col_k.push_back(coeff(r, r));
	}
	//TODO: if zero, need to halt linear search.
	
	if (coeff(k, k) !=0){
		col_r_nnzs.push_back(r);
		col_r.push_back(coeff(k, k));
	}
	
	for (i = 0; i < (int) list[r].size(); i++) {
		// if (list[r][i] < k) {  //list[r] does not update correctly on each iteration. the invariant is somehow violated. the code below works but is bad code. pretty sure this only fixes a subset of the problems.
			// std::swap(list[r][i], list[r][0]);
			// list[r].pop_front();
			// i--;
			// continue;
		// }
		
		if (list[r][i] < k) continue;
		coeffRef(r, list[r][i], its_k);
		col_k_nnzs.push_back(list[r][i]);
		col_k.push_back(*its_k.second);
		
		*its_k.first = m_idx[list[r][i]].back();
		*its_k.second = m_x[list[r][i]].back();
		
		m_idx[list[r][i]].pop_back();
		m_x[list[r][i]].pop_back();
	}
	
	if (coeff(r, k) != 0) {
		col_k_nnzs.push_back(r);
		col_k_nnzs.push_back(coeff(r, k));
	}


	//is it possible to select an all zeros col as a pivot?
	if (m_idx[r].size() > 0) {
		offset = (m_idx[r][0] == r ? 1 : 0);
		std::copy(m_x[r].begin()+offset, m_x[r].end(), std::back_inserter(col_k));
		std::copy(m_idx[r].begin()+offset, m_idx[r].end(), std::back_inserter(col_k_nnzs));
		
		//invariant: ensure list[:][0] contain the index nearest in value to k.
		for (auto it = m_idx[r].begin() + offset; it != m_idx[r].end(); it++) {
			for (i = 0; i < (int) list[*it].size(); i++) {
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
			if (m_idx[k][i] < r) {
				m_idx[m_idx[k][i]].push_back(r);
				m_x[m_idx[k][i]].push_back(m_x[k][i]);
			} else if (m_idx[k][i] != r) {
				col_r.push_back(m_x[k][i]);
				col_r_nnzs.push_back(m_idx[k][i]);
				
				//this part can be simplified since invariant ensures list[:][0] is always the idx closest to k. wont make a diff on running time though
				for (j = 0; j < (int) list[m_idx[k][i]].size(); j++) {
					if (list[m_idx[k][i]][j] == k) {
						swapr_.push_back(list[m_idx[k][i]].begin() + j);
						all_swaps.push_back(m_idx[k][i]);
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
		ensure_invariant(*it, k, list[*it]);
	}
	
	//set the kth col
	m_idx[k].swap(col_k_nnzs);
	m_x[k].swap(col_k);
	
	//set the rth col
	m_idx[r].swap(col_r_nnzs);
	m_x[r].swap(col_r);
	
	list[k].swap(list[r]);
	//--------end pivot A---------//
	
	//----------pivot L ----------//
	swapr.clear();
	swapk.clear();
	all_swaps.clear();
	
	for (auto it = L.list[k].begin(); it != L.list[k].end(); it++)
	{
		for (i = L.first[*it]; i < (int) L.m_idx[*it].size(); i++) {
			if (L.m_idx[*it][i] == k) {
				swapr.push_back(L.m_idx[*it].begin() + i);
				break;
			}
		}
	}
	
	for (auto it = L.list[r].begin(); it != L.list[r].end(); it++) {
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