// -*- mode: c++ -*-
#ifndef SQUARE_MATRIX_FIND_ROOT_H
#define SQUARE_MATRIX_FIND_ROOT_H

template <class el_type>
inline void square_matrix<el_type>::find_root(int& s, vector<int>& deg)
{
	vector<bool> visited(m_n_cols);
	vector<int> current_level, next_level;
	current_level.reserve(m_n_cols);
	next_level.reserve(m_n_cols);
	int ls_max = 0, ls, d, min_deg;
	
	while (true)
	{
		std::fill(visited.begin(), visited.end(), false);
		ls = 0;
		current_level.clear();
		visited[s] = true;
		current_level.push_back(s);
		
		while (true)
		{
			for (auto node_it = current_level.begin(); node_it != current_level.end(); node_it++)
			{
				for (auto it = list[*node_it].begin(); it != list[*node_it].end(); it++)
				{
					if (!visited[*it])
					{
						visited[*it] = true;
						next_level.push_back(*it);
					}
				}
				
				for (auto it = m_idx[*node_it].begin(); it != m_idx[*node_it].end(); it++)
				{
					if (!visited[*it])
					{
						visited[*it] = true;
						next_level.push_back(*it);
					}
				}
			}
			if (next_level.empty())
				break;
			ls++;
			next_level.swap(current_level);
			next_level.clear();
		}

		if (ls > ls_max)
		{
			ls_max = ls;
			min_deg = m_n_cols;
			for (auto it = current_level.begin(); it != current_level.end(); it++)
			{
				d = deg[*it];
				if (d < min_deg)
				{ //should consider tie breaking by index later if needed.
					min_deg = d;
					s = *it;
				}
			}
		}
		else
			break;
	}
}

#endif