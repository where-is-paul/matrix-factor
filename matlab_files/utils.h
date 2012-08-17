#ifndef _MEX_UTILS
#define _MEX_UTILS

#include "../source/lilc_matrix.h"
#include "../source/solver.h"
#include <mex.h>
#include <algorithm>
#include <queue>

class mex_utils
{	
	public:
		solver<double> solv;
		
		// validation methods
		bool is_double_scalar(const mxArray* a);
		bool is_vector(const mxArray* a);
		bool is_double_vector(const mxArray* a);
		bool is_string(const mxArray* a);

		// parsing methods
		double parse_double(const mxArray* raw_alpha);
		
		void mex_convert(double* m_x, mwSize* m_col_idx, mwSize* m_row_idx, mwSize& nnzs);
		
		void mex_set(double*& m_x, mwSize*& m_col_idx, mwSize*& m_row_idx, mwSize& nnzs, mwSize& m, mwSize& n, char matrix_type);
		
		void mex_save_lhs(mxArray*& lhs_ptr, char matrix_type);
		
		~mex_utils() {
		}		
};

#endif