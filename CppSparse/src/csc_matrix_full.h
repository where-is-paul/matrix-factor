//-*-mode:c++-*-
#ifndef _CSC_MATRIX_FULL_H_
#define _CSC_MATRIX_FULL_H_
template <class idx_type, class el_type>
std::vector<el_type> csc_matrix <idx_type, el_type> :: full () const
{
    std::vector<el_type> result (m_n_rows * m_n_cols);

    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        auto col_iter = result.begin() + m_n_rows * col;
        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {
            idx_type row = m_row_idx[offset];
            el_type  nzv = m_x [offset];
            *(col_iter + row) += nzv;
        }
    }
    return result;
}
#endif
