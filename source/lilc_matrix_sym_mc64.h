// -*- mode: c++ -*-
#ifndef _LILC_MATRIX_SYM_MC64_H_
#define _LILC_MATRIX_SYM_MC64_H_

#include <vector>

#include "hsl_mc64d.h"

// Change include file above to change support to floats, complex, complex doubles, etc.
template<> 
inline std::vector<double> lilc_matrix<double> :: sym_mc64(vector<int>& perm) {
	std::vector<double> row_val;
	std::vector<int> row_ind, col_ptr;
	
	to_csc(row_val, row_ind, col_ptr);
	
	const int matrix_type = 4; //Symmetric matrix type
	const int job = 5; // Symmetrized permutation and scaling. We dont use the scaling for now.
	int* ptr = &col_ptr[0];
	int* row = &row_ind[0];
	double* val = &row_val[0];
	
	struct mc64_control control;
	struct mc64_info info;
	
	// Initialize control
	mc64_default_control(&control);
	
	// Last null parameter is for the scaling, which we currently don't use.
	int* p = new int[m_n_rows + m_n_cols];
	double *s = new double[m_n_rows + m_n_cols];
	mc64_matching(job, matrix_type, m_n_rows, m_n_cols, ptr, row, val, &control, &info, p, s);
	if (info.flag < 0) {
		std::cerr << "Failure of mc64 with flag " << info.flag << std::endl;
	}
	
	perm.assign(p, p+m_n_cols);
	delete p;
	
	// Return scaling factors
	std::vector<double> res = std::vector<double>(s, s+m_n_cols);
	for (int i = 0; i < m_n_cols; i++) {
		res[i] = exp(res[i]);
	}
	delete s;
	
	return res;
}

#endif