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
    if (nlhs > 3)
        mexErrMsgTxt("Too many output arguments.");
    // set up raw variables
    const mxArray* raw_csc = prhs[0];
	const mxArray* raw_fill_fact = prhs[1];
    const mxArray* raw_tol = prhs[2];
	
	double fill_factor = parse_double(raw_fill_fact);
	double tol = parse_double(raw_tol);

    if (mxGetN(raw_csc) != mxGetM(raw_csc))
        mexErrMsgTxt("matrix must be square.");
	
	double* m_x, *L_x, *D_x;
	unsigned long int* m_row_idx, *L_row_idx, *D_row_idx;
	unsigned long int* m_col_idx, *L_col_idx, *D_col_idx;
	unsigned long int N, L_N, D_N;
	
	m_x = mxGetPr(raw_csc);
    m_col_idx = mxGetJc(raw_csc);
    m_row_idx = mxGetIr(raw_csc);
	N = mxGetN(raw_csc);
	
	solver<double> solv;
	solv.mex_convert(m_x, m_col_idx, m_row_idx, N);
	solv.solve(fill_factor, tol);
	
	//-------------- set L ---------------//
	plhs[0] =  mxCreateSparse((mwSize) solv.L.n_rows(), 
							  (mwSize) solv.L.n_cols(), 
							  (mwSize) solv.L.nnz(), 
							  (mxComplexity) NULL);
	L_x  = mxGetPr(plhs[0]);
    L_col_idx = mxGetJc(plhs[0]);
	L_row_idx = mxGetIr(plhs[0]);
	
	solv.mex_set_L(L_x, L_col_idx, L_row_idx, L_N);
	
	mxSetNzmax(plhs[0], L_N);
	mxSetPr(plhs[0], mxRealloc(L_x, L_N*sizeof(double)));
	mxSetIr(plhs[0], mxRealloc(L_row_idx, L_N*sizeof(mwIndex)));
	mxSetJc(plhs[0], mxRealloc(L_col_idx, L_N*sizeof(mwIndex)));
	//------------------------------------//
	
	//-------------- set D ---------------//
	plhs[1] =  mxCreateSparse((mwSize) solv.D.n_rows(), 
							  (mwSize) solv.D.n_cols(), 
							  (mwSize) solv.D.nnz(), 
							  (mxComplexity) NULL);
	D_x  = mxGetPr(plhs[1]);
    D_col_idx = mxGetJc(plhs[1]);
	D_row_idx = mxGetIr(plhs[1]);
	
	solv.mex_set_D(D_x, D_col_idx, D_row_idx, D_N);
	
	mxSetNzmax(plhs[1], D_N);
	mxSetPr(plhs[1], mxRealloc(L_x, D_N*sizeof(double)));
	mxSetIr(plhs[1], mxRealloc(L_row_idx, D_N*sizeof(mwIndex)));
	mxSetJc(plhs[1], mxRealloc(L_col_idx, D_N*sizeof(mwIndex)));
	//------------------------------------//
	
}