//-*- mode: c++ -*-
#ifndef SYMMETRY_MATRIX_PERMUTE_H
#define SYMMETRY_MATRIX_PERMUTE_H

template<class el_type>
inline void square_matrix<el_type>::permute(std::vector<int>& perm)
{
	vector<idx_vector_type> m_idx_new(m_n_cols);
	vector<elt_vector_type> m_x_new(m_n_cols);
	
	int i, j, pi, pj;
	el_type px;
	vector<int> pinv(m_n_cols);
	for (j = 0, i = nnz() / m_n_cols; j < m_n_cols; j++)
	{
		pinv[perm[j]] = j;
		list[j].clear();
		
		m_idx_new[j].reserve(i);
		m_x_new[j].reserve(i);
	}
	
	for (j = 0; j < m_n_cols; j++)
	{ //no need to use function call n_cols() every iter
		pj = pinv[j];
		
		for (i = 0; i < (int)m_idx[j].size(); i++)
		{
			pi = pinv[m_idx[j][i]];
			px = m_x[j][i];
			
			if (pi < pj)
			{
				m_idx_new[pi].push_back(pj);
				m_x_new[pi].push_back(px * sign);
				list[pj].push_back(pi);
			}
			else
			{
				m_idx_new[pj].push_back(pi);
				m_x_new[pj].push_back(px);
				if (pi != pj) list[pi].push_back(pj);
			}
		}
	}
	//vector<int> count(50, 0);
	//for (int i = 0; i < m_n_cols; i++)
	//	count[m_idx_new[i].size()]++;
	//for (int i = 0; i < (int)count.size(); i++)
	//	printf("Size %d: %d.\n", i, count[i]);
	
	m_idx.swap(m_idx_new);
	m_x.swap(m_x_new);
}

#endif