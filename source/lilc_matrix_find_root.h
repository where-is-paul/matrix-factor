// -*- mode: c++ -*-
#ifndef _LILC_MATRIX_FIND_ROOT_H_
#define _LILC_MATRIX_FIND_ROOT_H_

template<class el_type> 
void lilc_matrix<el_type> :: find_root(int& s) {
	vector<bool> visited(m_n_cols, false);
	vector<int> lvl_set;
	int ls = 0;
	
	lvl_set.push_back(s);
	visited[s] = true;
	while (find_level_set(lvl_set, visited))
		ls++;
	
	int deg, min_deg = m_n_cols;
	for (auto it = lvl_set.begin(); it != lvl_set.end(); it++) {
		deg = list[*it].size() + m_idx[*it].size();
		if (deg < min_deg) { //should consider tie breaking by index later if needed.
			min_deg = deg;
			s = *it;
		}
	}
}

#endif