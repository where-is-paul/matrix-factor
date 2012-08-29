#include "utils.h"

using namespace std;

typedef pair<mwSize, double> idx_val_pair;

void load_index_value_pairs(vector<idx_val_pair>& v12, vector<int>& v1, vector<double>& v2) {
	v12.clear();
	for (int i = 0; i < (int) v1.size(); i++) {
		v12.push_back(make_pair(v1[i], v2[i]));
	}
}

bool mex_utils :: is_double_scalar(const mxArray* a) {
    return mxIsDouble(a)
            && !mxIsComplex(a)
            && mxGetNumberOfElements(a) == 1;
}

bool mex_utils :: is_vector(const mxArray* a) {
    const mwSize* dims = mxGetDimensions(a);
    return (mxGetNumberOfDimensions(a) == 2
            && min(dims[0], dims[1]) <= 1);
}

bool mex_utils :: is_double_vector(const mxArray* a) {
    return mxIsDouble(a)
            && !mxIsComplex(a)
            && is_vector(a);
}

bool mex_utils :: is_string(const mxArray* a) {
    return mxIsChar(a) && is_vector(a);
}

double mex_utils :: parse_double(const mxArray* raw_double) {
    if (!is_double_scalar(raw_double))
        mexErrMsgTxt("Argument must be a real scalar.");
    double res = *((double*) mxGetData(raw_double));

    return res;
}

/*! \brief Loads the matrix A into solver (A is stored in CSC form).
	\param m_x the non-zero values stored in the matrix.
	\param m_col_idx the column pointers for the CSC matrix.
	\param m_row_idx the row pointers for the CSC matrix.
*/
void mex_utils :: mex_convert(double* m_x, mwSize* m_col_idx, mwSize* m_row_idx, mwSize& nnzs) {
	int count = 0;

	solv.A.resize(nnzs, nnzs);
	fill(solv.A.first.begin(), solv.A.first.end(), 0);
		
	for (mwSize i = 0; i < nnzs; i++) {

		for (mwSize j = m_col_idx[i]; j < m_col_idx[i+1]; j++) {

			if (m_row_idx[j] < i) continue;
			solv.A.m_idx[i].push_back(m_row_idx[j]);
			solv.A.m_x[i].push_back(m_x[j]);
			if (i != m_row_idx[j]) 
				solv.A.list[ m_row_idx[j] ].push_back(i);
			count++;
		}
	}
	
	solv.A.nnz_count = count;
	

}

void mex_utils :: mex_set(double*& m_x, mwSize*& m_col_idx, mwSize*& m_row_idx, mwSize& nnzs, mwSize& m, mwSize& n, char matrix_type) {
	m = solv.A.n_rows();
	n = solv.A.n_cols();
	int i = 0, count = 0;
	vector<idx_val_pair> temp_col;
	vector<queue<idx_val_pair> > elems(solv.A.n_cols());
	idx_val_pair elem;
	
	switch (matrix_type) {
		case 'A':
			m_x = (double*) mxCalloc(2*solv.A.nnz(), sizeof(double));
			m_row_idx = (mwSize*) mxCalloc(2*solv.A.nnz(), sizeof(mwSize));
			m_col_idx = (mwSize*) mxCalloc(solv.A.n_cols()+1, sizeof(mwSize));
			
			for (i = 0; i < (int) solv.A.n_rows(); i++) {
				m_col_idx[i] = count;
				while (!elems[i].empty() ) {
					elem = elems[i].front();
					elems[i].pop();
					
					m_row_idx[count] = elem.first;
					m_x[count] = elem.second;
					count++;
				}
				
				load_index_value_pairs(temp_col, solv.A.m_idx[i], solv.A.m_x[i]);
				sort(temp_col.begin(), temp_col.end());
				
				for (unsigned int j = 0; j < temp_col.size(); j++) {
					m_row_idx[count] = temp_col[j].first;
					m_x[count] = temp_col[j].second;
					count++;
					
					if (i != (int) temp_col[j].first)
					elems[ temp_col[j].first ].push( 
						make_pair(i, temp_col[j].second)
					);
				}
			}
			
			break;
			
		case 'L':
			m_x = (double*) mxCalloc(solv.L.nnz(), sizeof(double));
			m_row_idx = (mwSize*) mxCalloc(solv.L.nnz(), sizeof(mwSize));
			m_col_idx = (mwSize*) mxCalloc(solv.L.n_cols()+1, sizeof(mwSize));
			for (i = 0; i < (int) solv.L.n_rows(); i++) {
				m_col_idx[i] = count;
				
				load_index_value_pairs(temp_col, solv.L.m_idx[i], solv.L.m_x[i]);
				sort(temp_col.begin(), temp_col.end());
				
				for (unsigned int j = 0; j < temp_col.size(); j++) {
					m_row_idx[count] = temp_col[j].first;
					m_x[count] = temp_col[j].second;
					count++;
				}
			}
			
			break;
		
		case 'D':
			m_x = (double*) mxCalloc(2*solv.D.nnz(), sizeof(double));
			m_row_idx = (mwSize*) mxCalloc(2*solv.D.nnz(), sizeof(mwSize));
			m_col_idx = (mwSize*) mxCalloc(solv.D.n_cols()+1, sizeof(mwSize));
			
			for (i = 0; i < solv.D.n_rows(); i++) {
				m_col_idx[i] = count;
				while (!elems[i].empty() ) {
					elem = elems[i].front();
					elems[i].pop();
					
					m_row_idx[count] = elem.first;
					m_x[count] = elem.second;
					count++;
				}
				
				m_row_idx[count] = i;
				m_x[count] = solv.D[i];
				if (solv.D.block_size(i) == 2) {
					count++;
					m_row_idx[count] = i+1;
					m_x[count] = solv.D.off_diagonal(i);
					
					elems[ i+1 ].push( 
						make_pair(i, solv.D.off_diagonal(i)) 
					);
					
					
				}
				count++;
			}
			break;
		
		case 'S':
			m_x = (double*) mxCalloc(solv.A.S.nnz(), sizeof(double));
			m_row_idx = (mwSize*) mxCalloc(solv.A.S.nnz(), sizeof(mwSize));
			m_col_idx = (mwSize*) mxCalloc(solv.A.S.n_cols()+1, sizeof(mwSize));
			for (i = 0; i < solv.A.S.n_rows(); i++) {
				m_col_idx[i] = count;
				m_row_idx[count] = i;
				m_x[count] = solv.A.S[i];
				count++;
			}
			break;
		
		case 'P':
			n = 1;
			m_x = (double*) mxCalloc(solv.perm.size(), sizeof(double));
			m_row_idx = (mwSize*) mxCalloc(solv.perm.size(), sizeof(mwSize));
			m_col_idx = (mwSize*) mxCalloc(2, sizeof(mwSize));
			
			m_col_idx[i] = count;
			for (int j = 0; j < (int) solv.perm.size(); j++) {
				m_row_idx[count] = j;
				m_x[count] = solv.perm[j]+1;
				count++;
			}
			i++;
			
			break;
		
		default:
			cerr << "No options selected for mex_set." << endl;
	}
	
	m_col_idx[i] = count;
	nnzs = count;
}

void mex_utils :: mex_save_lhs(mxArray*& lhs_ptr, char matrix_type) {
	double* m_x;
	mwSize* m_row_idx;
	mwSize* m_col_idx;
	mwSize nnzs, m, n;
	
	lhs_ptr =  mxCreateSparse(0, 0, 0, mxREAL);
							  
	m_x  = mxGetPr(lhs_ptr);
	m_row_idx = mxGetIr(lhs_ptr);
    m_col_idx = mxGetJc(lhs_ptr);
	
	mex_set(m_x, m_col_idx, m_row_idx, nnzs, m, n, matrix_type);
	
	mxSetM(lhs_ptr, m);
	mxSetN(lhs_ptr, n);
	mxSetNzmax(lhs_ptr, nnzs);
	mxSetPr(lhs_ptr, m_x);
	mxSetIr(lhs_ptr, m_row_idx);
	mxSetJc(lhs_ptr, m_col_idx);
	mxSetPi(lhs_ptr, NULL);
}