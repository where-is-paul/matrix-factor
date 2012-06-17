// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_REACH_H_
#define _CSC_MATRIX_REACH_H_

#include "flip.h"

template<class idx_type, class el_type>
idx_type csc_matrix<idx_type, el_type>::reach(const csc_t& b, idx_type col, const idx_vector_type& pinv, idx_vector_type& xi) const
{
    idx_type top = m_n_cols;
    idx_vector_type pstack(m_n_cols);
    for (idx_type offset = b.m_col_idx[col]; offset < b.m_col_idx[col + 1]; offset ++)
    {
        idx_type row = b.m_row_idx[offset];

        if (!(is_flipped(m_col_idx[row])))
        {
            dfs(row, xi, pstack, pinv, top);
        }
    }

    auto col_idx = const_cast<idx_vector_type&>(m_col_idx);

    for (idx_type offset = top; offset < m_n_cols; offset ++)
    {
        col_idx[xi[offset]] = flip(col_idx[xi[offset]]);
    }

    return top;
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type>::dfs(idx_type col, idx_vector_type& xi, 
                                        idx_vector_type& pstack, const idx_vector_type& pinv, idx_type& top) const
{
    idx_type head = 0;
    xi[head]      = col;
    
    while (head != npos)
    {
        idx_type col = xi[head];

        if (!is_flipped(m_col_idx[col]))
        {
            pstack[head]   = m_col_idx[col];
            idx_vector_type& col_idx = const_cast<idx_vector_type&>(m_col_idx);
            col_idx[col] = flip(col_idx[col]);
        }

        bool done      = true;
        idx_type pstop = unflip(m_col_idx[col+1]);
        for (idx_type offset = pstack[head]; offset < pstop; offset ++)
        {
            idx_type row =  m_row_idx[offset];
            if (!is_flipped(m_col_idx[row]))
            {
                pstack[head] = offset;
                xi[++head]   = row;
                done         = false;
                break;
            }
        }

        if (done)
        {
            --head;
            xi[--top] = col;
        }
    }
} 

#endif
