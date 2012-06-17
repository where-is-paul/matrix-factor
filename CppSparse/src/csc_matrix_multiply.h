// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_MULTIPLY_H_
#define _CSC_MATRIX_MULTIPLY_H_


template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: multiply (const csc_matrix<idx_type, el_type>& b) const
{
    csc_t result (n_rows(), b.n_cols(), nnz() + b.nnz()); 

    idx_type        nnz_result = 0;
    idx_vector_type work (result.n_rows());
    elt_vector_type x    (result.n_rows());

    for (idx_type j = 0; j < result.m_n_cols; j ++)
    {
        for (idx_type offset = b.m_col_idx[j]; offset < b.m_col_idx[j+1]; offset ++)
        {
            idx_type k    = b.m_row_idx[offset];
            el_type  beta = b.m_x[offset];
            scatter (result, k, beta, work, x, j+1, nnz_result);
        }
        result.m_col_idx[j+1] = nnz_result;
        result.gather (j, x);
    }
    result.m_row_idx.resize(nnz_result);
    result.m_x.resize(nnz_result);
    return result;
}



#endif // _CSC_MATRIX_MULTIPLY_H_
