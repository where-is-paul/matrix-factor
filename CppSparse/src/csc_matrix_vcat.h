//-*- mode: c++ -*-
#ifndef _CSC_MATRIX_VCAT_H_
#define _CSC_MATRIX_VCAT_H_


template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: do_vcat (idx_type row_offset, std::vector<idx_type>& work, 
                                               const csc_matrix<idx_type, el_type>& x)
{
    const std::vector<idx_type>& x_col_idx = x.col_idx();
    const std::vector<idx_type>& x_row_idx = x.row_idx();
    const std::vector<el_type >& x_data    = x.nz_vals();

    for (idx_type col = 0; col < x.m_n_cols; col ++)
    {
        for (idx_type offset = x_col_idx[col]; offset < x_col_idx[col + 1]; offset ++)
        {
            idx_type& real_offset  = work[col];
            m_row_idx[real_offset] = x_row_idx[offset] + row_offset;
            m_x      [real_offset] = x_data [offset];
            real_offset ++;
        }
    }
}


#endif
