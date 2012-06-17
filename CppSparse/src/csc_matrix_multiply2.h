// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_MULTIPLY2_H_
#define _CSC_MATRIX_MULTIPLY2_H_


template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: multiply2 (const csc_matrix<idx_type, el_type>& b) const
{
    idx_type m = n_rows();
    idx_type k = n_cols();
    idx_type n = b.n_cols();

    csc_t c(m, n);

    idx_vector_type new_offsets(m, npos);

    idx_type new_offset = 0;

    for (idx_type j = 0; j < n; j ++)
    {
        idx_type col_start = new_offset;

        for (idx_type offset_b = b.m_col_idx[j]; offset_b < b.m_col_idx[j + 1]; offset_b ++)
        {
            idx_type k = b.m_row_idx[offset_b];
            for (idx_type offset_a =  m_col_idx[k]; offset_a < m_col_idx[k + 1]; offset_a ++)
            {
                idx_type i = m_row_idx[offset_a];
                if (new_offsets[i] == npos || new_offsets[i] < col_start)
                {
                    new_offsets[i] = new_offset;
                    new_offset ++;
                }
            }
        }
        
        c.m_col_idx[j] = col_start;
    }
    
    c.m_col_idx.back() = new_offset;    
    c.m_row_idx.resize(new_offset);
    c.m_x.resize(new_offset);
    
    new_offset = 0; 
    std:fill(new_offsets.begin(), new_offsets.end(), npos);

    for (idx_type j = 0; j < n; j ++)
    {
        idx_type col_start = new_offset;

        for (idx_type offset_b = b.m_col_idx[j]; offset_b < b.m_col_idx[j + 1]; offset_b ++)
        {
            idx_type k  = b.m_row_idx[offset_b];
            el_type  bx = b.m_x      [offset_b];
            for (idx_type offset_a =  m_col_idx[k]; offset_a < m_col_idx[k + 1]; offset_a ++)
            {
                idx_type i  = m_row_idx[offset_a];
                el_type  ax = m_x      [offset_a];
                if (new_offsets[i] == npos || new_offsets[i] < col_start)
                {
                    new_offsets[i] = new_offset;
                    c.m_row_idx[new_offset] = i;
                    c.m_x      [new_offset] = ax*bx;
                    new_offset ++;
                }
                else
                {
                    c.m_x [new_offsets[i]] += ax * bx;
                }
            }
        }
    }

    return c;
}



#endif // _CSC_MATRIX_MULTIPLY2_H_
