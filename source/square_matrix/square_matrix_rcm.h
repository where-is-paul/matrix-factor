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
	vector<int> q(m_n_cols+1);
	int head = 0, tail = 0;
	
	// precalculate deg and adjacency_list because they are used a large number of times
	vector<int> deg(m_n_cols);
	for (i = 0; i < m_n_cols; i++)
		deg[i] = degree(i);
	
	auto sorter = [&deg](const int &a, const int &b) -> bool
						{
							if (deg[a] == deg[b]) return a > b;
							return deg[a] < deg[b];
						};
						
	vector<int> nbrs;
	nbrs.reserve(m_n_cols);
	
	vector<int> index(m_n_cols);
	for (i = 0; i < m_n_cols; i++) index[i] = i;
	
	//sort(index.begin(), index.end(), sorter);
	//std::random_shuffle(index.begin(), index.end());
	//swap(index[22965], index[0]);
	for (i = 0; i < m_n_cols; i++)
	{
		if (visited[i])
			continue;
		
		s = index[i];
		find_root(s, deg);
		visited[s] = true;
		q[tail] = s;
		tail++;
		
		while (head != tail)
		{
			s = q[head];
			head++;
			perm.push_back(s);
			
			nbrs.clear();
			for (auto it = list[s].begin(); it != list[s].end(); it++)
			{
				if (!visited[*it]) {
					visited[*it] = true;
					nbrs.push_back(*it);
				}
			}
			
			for (auto it = m_idx[s].begin(); it != m_idx[s].end(); it++) {
				if (!visited[*it]) {
					visited[*it] = true;
					nbrs.push_back(*it);
				}
			}
			
			sort(nbrs.begin(), nbrs.end(), sorter);
			
			for (auto it = nbrs.begin(); it != nbrs.end(); it++) {
				q[tail++] = *it;
			}
		}
	}
	
	reverse(perm.begin(), perm.end());
}

#endif