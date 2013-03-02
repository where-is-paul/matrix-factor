// -*- mode: c++ -*-
#ifndef SQUARE_MATRIX_RCM_H
#define SQUARE_MATRIX_RCM_H

//#include <queue>
//using std::queue;

template<class el_type>
inline void square_matrix<el_type>::rcm(vector<int>& perm)
{
	int i, s, head = 0, tail = 0;
	vector<bool> visited(m_n_cols, false);
	//queue<int> q;
	vector<int> q(m_n_cols);
	vector<int> neighbors;
	neighbors.reserve(m_n_cols);

	// precalculate deg and adjacency_list because they are used a large number of times
	vector<int> deg(m_n_cols);
	for (i = 0; i < m_n_cols; i++)
		deg[i] = degree(i);
	
	auto sorter = [&deg](const int &a, const int &b) -> bool
						{
							if (deg[a] == deg[b]) return a > b;
							return deg[a] < deg[b];
						};

	for (i = 0; i < m_n_cols; i++)
	{
		if (visited[i])
			continue;
		
		s = i;
		find_root(s, deg);
		visited[s] = true;
		q[tail] = s; tail++; // equivalent to q.push(s);
		
		while (head != tail) // equivalent to while (!q.empty())
		{
			s = q[head]; head++;// equivalent to s = q.front(); q.pop();
			perm.push_back(s);
			neighbors.clear();

			for (auto it = list[s].begin(); it != list[s].end(); it++)
			{
				if (!visited[*it])
				{
					visited[*it] = true;
					neighbors.push_back(*it);
				}
			}
			for (auto it = m_idx[s].begin(); it != m_idx[s].end(); it++)
			{
				if (!visited[*it])
				{
					visited[*it] = true;
					neighbors.push_back(*it);
				}
			}
			sort(neighbors.begin(), neighbors.end(), sorter);

			for (auto it = neighbors.begin(); it != neighbors.end(); it++)
			{
				q[tail] = *it; tail++; // equivalent to q.push(*it);
			}
		}
	}
	
	reverse(perm.begin(), perm.end());
}

#endif