//-*-mode:c++-*-
#ifndef _CSC_MATRIX_FIND_H_
#define _CSC_MATRIX_FIND_H_

template<class idx_type, class el_type>
triplet_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: find () const
{
    triplet_matrix<idx_type, el_type> result (m_n_rows, m_n_cols, nnz());
    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        for (idx_type offset = m_col_idx [col]; offset < m_col_idx [col + 1]; offset ++)
        {
            result.push_back(m_row_idx[offset], col, m_x[offset]);
        }
    }
    
    return result;
}

#endif
