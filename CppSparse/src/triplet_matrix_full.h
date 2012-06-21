//-*-mode:c++-*-
#ifndef _TRIPLET_MATRIX_FULL_H_
#define _TRIPLET_MATRIX_FULL_H_

template <class idx_type, class el_type>
std::vector<el_type> triplet_matrix<idx_type, el_type> :: full () const
{
    std::vector<el_type> result (m_n_rows * m_n_cols);
    std::fill (result.begin(), result.end(), static_cast<el_type>(0));
    for (idx_type elem = 0; elem < m_x.size(); elem ++)
    {
        idx_type row = m_row_idx[elem];
        idx_type col = m_col_idx[elem];
        el_type  nzv = m_x [elem];
        result[row + m_n_rows * col] += nzv;
    }
    return result;
}

#endif
