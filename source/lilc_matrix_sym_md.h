// -*- mode: c++ -*-
#ifndef _LILC_MATRIX_SYM_MD_H_
#define _LILC_MATRIX_SYM_MD_H_

#include "quotient_graph.h"

template<class el_type> 
inline void lilc_matrix<el_type> :: sym_md(vector<int>& perm) {
	quotient_graph g(m_idx);
	
	int min_node; perm.resize(m_n_cols);
	for (int i = 0; i < m_n_cols; i++) {
		min_node = g.min_deg();
		g.eliminate(min_node);
		perm[i] = min_node;
	}
	
}

#endif