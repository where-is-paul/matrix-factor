//-*-mode:c++-*-
#ifndef _CSC_MATRIX_SCALE_H_
#define _CSC_MATRIX_SCALE_H_
#include <stdexcept>

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type>::scale (const std::vector<el_type>& r, const std::vector<el_type>& c) 
{
    if (r.size() < n_rows() || c.size() < n_cols())
    {
        throw std::out_of_range("scale");
    }

    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col+1]; offset ++)
        {
            m_x[offset] = r[m_row_idx[offset]] * m_x[offset] * c[col];
        }
    }
}
#endif
