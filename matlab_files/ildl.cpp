#include <mex.h>
#include <string>
#include "../source/solver.h"
#include "utils.h"

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // validate number of inputs and outputs
    // validate number of inputs and outputs
    if (nrhs < 1)
        mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 3)
        mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 5)
        mexErrMsgTxt("Too many output arguments.");
		
	mex_utils m;
	
    // set up raw variables
    const mxArray* raw_csc = prhs[0];
	
	double fill_factor = 1.0, tol = 0.001;
	if (nrhs > 1) {
		const mxArray* raw_fill_fact = prhs[1];
		fill_factor = m.parse_double(raw_fill_fact);
	}
	if (nrhs > 2) {
		const mxArray* raw_tol = prhs[2];
		tol = m.parse_double(raw_tol);
	}

    if (mxGetN(raw_csc) != mxGetM(raw_csc))
        mexErrMsgTxt("matrix must be square.");
	
	//--------------- load A ---------------//
	double* m_x;
	mwSize* m_row_idx;
	mwSize* m_col_idx;
	mwSize nnzs;
	
	m_x = mxGetPr(raw_csc);
    m_col_idx = mxGetJc(raw_csc);
    m_row_idx = mxGetIr(raw_csc);
	nnzs = mxGetN(raw_csc);
	
	m.mex_convert(m_x, m_col_idx, m_row_idx, nnzs);
	//--------------------------------------//
	
	//factor A.
	m.solv.solve(fill_factor, tol);

	m.mex_save_lhs(plhs[0], 'L');
	m.mex_save_lhs(plhs[1], 'D');
	m.mex_save_lhs(plhs[2], 'P');
	m.mex_save_lhs(plhs[3], 'S');
	
	if (nlhs == 5)
	m.mex_save_lhs(plhs[4], 'A');

}