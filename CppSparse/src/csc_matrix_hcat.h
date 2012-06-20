//-*- mode: c++ -*-
#ifndef _CSC_MATRIX_HCAT_H_
#define _CSC_MATRIX_HCAT_H_


template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: do_hcat (idx_type col_offset, const csc_matrix<idx_type, el_type>& x)
{
    const std::vector<idx_type>& x_col_idx = x.col_idx();
    const std::vector<idx_type>& x_row_idx = x.row_idx();
    const std::vector<el_type >& x_data    = x.nz_vals();

    const idx_type cum_offset = m_col_idx[col_offset];

    for (idx_type col = 1; col < x.n_cols() + 1; col ++)
    {
        m_col_idx [col_offset + col] = x_col_idx[col] + cum_offset;
    }

    std::copy(x_row_idx.cbegin(), x_row_idx.cend(), m_row_idx.begin() + m_col_idx[col_offset]);
    std::copy(x_data.cbegin()   , x_data.cend()   , m_x.begin()       + m_col_idx[col_offset]);
}

#endif
