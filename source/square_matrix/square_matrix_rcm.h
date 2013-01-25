// -*- mode: c++ -*-
#ifndef SQUARE_MATRIX_RCM_H
#define SQUARE_MATRIX_RCM_H

#include <queue>

using std::queue;

template<class el_type>
inline void square_matrix<el_type>::rcm(vector<int>& perm)
{
	int i, s;
	vector<bool> visited(m_n_cols, false);
	queue<int> q;
	
	// precalculate deg and adjacency_list because they are used a large number of times
	vector<int> deg(m_n_cols);
	for (i = 0; i < m_n_cols; i++)
		deg[i] = degree(i);
	
	auto sorter = [&deg](const int &a, const int &b) -> bool
						{
							if (deg[a] == deg[b]) return a > b;
							return deg[a] < deg[b];
						};
						
	vector<idx_vector_type> adjacency_list(m_n_cols);
	for (i = 0; i < m_n_cols; i++)
	{
		adjacency_list[i].reserve(deg[i]);
		adjacency_list[i].assign(list[i].begin(), list[i].end());
		if (!m_idx[i].empty())
		{
			s = (m_idx[i][0] == i ? 1 : 0);
			std::copy(m_idx[i].begin() + s, m_idx[i].end(), std::back_inserter(adjacency_list[i]));
		}
		sort(adjacency_list[i].begin(), adjacency_list[i].end(), sorter);
	}

	for (i = 0; i < m_n_cols; i++)
	{
		if (visited[i])
			continue;
		
		s = i;
		find_root(s, adjacency_list, deg);
		visited[s] = true;
		q.push(s);
		
		while (!q.empty())
		{
			s = q.front();
			q.pop();
			perm.push_back(s);
			for (auto it = adjacency_list[s].begin(); it != adjacency_list[s].end(); it++)
			{
				if (!visited[*it])
				{
					visited[*it] = true;
					q.push(*it);
				}
			}
		}
	}
	
	reverse(perm.begin(), perm.end());
}

#endif