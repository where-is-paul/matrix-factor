// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_TRIANGULAR_SOLVE_H_
#define _CSC_MATRIX_TRIANGULAR_SOLVE_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: lsolve  (elt_vector_type& b) const
{
    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        idx_type offset = m_col_idx[col];
        assert(m_row_idx[offset] == col);

        /**
         *   x(j) /= A(j, j)
         */
        b[col] /= m_x[offset];
        offset ++;
        el_type ajj = b[col];
        for (offset; offset < m_col_idx[col + 1]; offset ++)
        {
            /**
             *  for i = j + 1 : n
             *     x(i) -= A(i, j) * x(j)
             *  end
             */
            b[m_row_idx[offset]] -= m_x[offset] * ajj;
        }
    }
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: ltsolve (elt_vector_type& b) const
{
    idx_type col = 0, j = m_n_cols - 1;
    for (; col < m_n_cols; col ++, j --)
    {
        idx_type offset = m_col_idx[j];
        idx_type row    = m_row_idx[offset];
        el_type  ajj    = m_x[offset];
        assert(row == j);

        /**
         * for i = j + 1 : n
         *    b(j) -= A(i, j) * b(j)
         * end
         */
        for (offset ++; offset < m_col_idx[j + 1]; offset ++)
        {
            idx_type row = m_row_idx[offset];
            el_type  aij = m_x[offset];
            b[j] -= aij * b[row];
        }

        /**
         * b(j) /= A(j, j)
         */
        b[j] /= ajj;
    }
}


template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: usolve  (elt_vector_type& b) const
{
    idx_type col = 0, j = m_n_cols - 1;
    for (idx_type col = 0; col < m_n_cols; col ++, j--)
    {
        idx_type diag_offset = m_col_idx[j + 1] - 1;        
        assert(m_row_idx[diag_offset] == j);
        b[j] /= m_x[diag_offset];
        for (idx_type offset = m_col_idx[j]; offset < diag_offset; offset ++)
        {
            idx_type row = m_row_idx[offset];
            el_type  aij = m_x[offset];
            b[row] -= aij * b[j];
        }
    }
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: utsolve (elt_vector_type& b) const
{
    for (idx_type j  = 0; j < m_n_cols; j ++)
    {
        idx_type offset = m_col_idx[j];
        for (; offset < m_col_idx[j + 1] - 1; offset ++)
        {
            el_type  aij = m_x      [offset];
            idx_type i   = m_row_idx[offset];
            b[j]        -= aij * b[i];
        }
        /**
         * Last element in the row indices is the diagonal element
         */
        assert(m_row_idx[offset] == j);
        b[j] /= m_x[offset];
    }
}



#endif  // _CSC_MATRIX_TRIANGULAR_SOLVE_H_
