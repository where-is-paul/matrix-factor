// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_ADD_H_
#define _CSC_MATRIX_ADD_H_


template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: add (const csc_matrix<idx_type, el_type>& b, 
                                                                    el_type alpha, el_type beta) const
{
    if (shape() != b.shape())
        throw std::logic_error("Shape mismatch in csc_matrix_t::add");

    // For A + B nnz(A) + nnz(B) is an overestimate for the number of
    // non-zero entries
    csc_t result(m_n_rows, m_n_cols, nnz() + b.nnz());
    
    idx_type nnz_result = 0;

    const csc_t& a = *this;

    idx_vector_type work (m_n_rows);
    elt_vector_type x    (m_n_rows);
    
    for (idx_type j = 0; j < m_n_cols; j ++)
    {
        a.scatter (result, j, alpha, work, x, j+1, nnz_result);
        b.scatter (result, j, beta , work, x, j+1, nnz_result);
        result.m_col_idx[j+1] = nnz_result;
        result.gather (j, x);
    }

    result.m_row_idx.resize(nnz_result);       
    return result;
}

#endif
