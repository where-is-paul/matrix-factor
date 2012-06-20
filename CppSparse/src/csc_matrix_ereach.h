// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_EREACH_H_
#define _CSC_MATRIX_EREACH_H_
#include "flip.h"

template<class idx_type, class el_type>
idx_type csc_matrix<idx_type, el_type>::ereach(const idx_vector_type& parent, 
                                               idx_type col, idx_vector_type& stack, idx_vector_type& work) const
{
    work[col] = flip(work[col]);
    idx_type top = m_n_rows;
    for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
    {
        idx_type row  = m_row_idx[offset];

        // Row indices of A are sorted
        if (row > col)
        {
            break;
        }

        // head contains the head of the queue used to traverse the
        // parent linked list Note that we need the reach set to be in
        // LIFO order in order to be considered "toplogically sorted"
        // but traversing the parent array returns the set in FIFO
        // order, so there in an additional reverse step at the end

        idx_type head = 0;

        // NOTE: We do not have to check for parent[row] == -1
        // NOTE: since by the very fact that row is in L[k,i] (since row is in A[k,i])
        // NOTE: at least k is guaranteed to be a parent of row and therefore
        // NOTE: cannot be -1
        while(!is_flipped(work[row]))
        {
            stack[head++] = row;
            work[row] = flip(work[row]);
            row = parent[row];
        }

        for (idx_type elem = 0; elem < head; elem ++)
        {
            stack[--top] = stack[elem];
        }
    }

    // Unflip all the elements of work
    for (idx_type elem = top; elem != m_n_rows; elem ++)
    {
        work[stack[elem]] = unflip(work[stack[elem]]);
    }

    return top;
}

template<class idx_type, class el_type>
typename csc_matrix<idx_type, el_type>::idx_vector_type csc_matrix<idx_type, el_type>::ereach(idx_type col) const
{
    if (m_n_rows != m_n_cols)
    {
        throw std::logic_error("Matrix must be square");
    }

    idx_vector_type parent = etree();

    idx_vector_type tstack(m_n_rows);
    idx_vector_type work  (m_n_rows);

    idx_type top = ereach(parent, col, tstack, work);
    idx_vector_type reach(m_n_rows - top);
    std::copy(tstack.cend() - top, tstack.cend(), reach.begin());
    return reach;
}

#endif // _CSC_MATRIX_EREACH_H_
