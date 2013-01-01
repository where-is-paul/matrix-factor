// -*- mode: c++ -*-
#ifndef HALF_MATRIX_RCM_H
#define HALF_MATRIX_RCM_H

template<class el_type> 
inline void half_matrix<el_type>::rcm(vector<int>& perm)
{
	int i, s;
	vector<bool> visited(m_n_cols, false);
	vector<int> lvl_set;
	for (i = 0; i < m_n_cols; i++)
	{
		if (visited[i])
			continue;
		
		lvl_set.clear();
		s = i;
		find_root(s);
		lvl_set.push_back(s);
		perm.push_back(s);
		visited[s] = true;

		auto sorter = [this](int const &a, int const &b) -> bool
						{
							int deg_a = this->degree(a), deg_b = this->degree(b);
							if (deg_a == deg_b) return a > b;
							return deg_a < deg_b;
						};
		
		while (find_level_set(lvl_set, visited))
		{
			sort(lvl_set.begin(), lvl_set.end(), sorter);
			perm.insert(perm.end(), lvl_set.begin(), lvl_set.end());
		}
	}
	
	reverse(perm.begin(), perm.end());
}

#endif