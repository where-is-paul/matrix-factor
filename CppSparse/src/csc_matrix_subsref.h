// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_SUBSREF_H_
#define _CSC_MATRIX_SUBSREF_H_

#ifdef __GNUG__
#define __FUNCDNAME__ __func__
#endif

template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> 
csc_matrix<idx_type, el_type>::subsref_range(idx_type i1, idx_type i2, idx_type j1, idx_type j2) const
{
    if (i1 >= m_n_rows || i2 >= m_n_rows || j1 >= m_n_cols || j2 >= m_n_cols)
    {
        throw std::out_of_range(__FUNCDNAME__);
    }

    csc_t result(i2-i1, j2-j1);
    idx_type curr_nnz = 0;
    for (idx_type col = j1; col < j2; col ++)
    {
        result.m_col_idx[col-j1] = curr_nnz;

        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {
            idx_type row = m_row_idx[offset];
            if (row >= i2)
            {
                break;
            }
            if (row >= i1)
            {
                curr_nnz ++;
            }
        }
    }
    result.m_col_idx.back() = curr_nnz;
    result.m_row_idx.resize(curr_nnz);
    result.m_x.resize(curr_nnz);
    
    curr_nnz = 0;
    for (idx_type col = j1; col < j2; col ++)
    {
        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {
            idx_type row = m_row_idx[offset];
            if (row >= i2)
            {
                break;
            }
            if (row >= i1)
            {
                result.m_row_idx[curr_nnz] = row;
                result.m_x      [curr_nnz] = m_x [offset];
                curr_nnz ++;
            }
        }
    }

    return result;
}


template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> 
csc_matrix<idx_type, el_type>::subsref_vector (const idx_vector_type& rows, const idx_vector_type& cols) const

{
    idx_type nr = n_rows();
    idx_type nc = n_cols();

    if (std::any_of(rows.cbegin(), rows.cend(), [nr](idx_type row) {return row >= nr;}))
    {
        throw std::out_of_range(__FUNCDNAME__);
    }

    if (std::any_of(cols.cbegin(), cols.cend(), [nc](idx_type col) {return col >= nc;}))
    {
        throw std::out_of_range(__FUNCDNAME__);
    }
    
    csc_t result(rows.size(), cols.size());
    
    /**
     * Subsref with random subscripts requires a binary search through
     * the row indices. Instead of performing a symbolic analysis to
     * determine the non-zero pattern, the better option is to just
     * append the new row indices and rely on the asymptotic O(1)
     * insertion into std::vector
     */

    for (auto col_iter = cols.cbegin(); col_iter < cols.cend(); col_iter ++)
    {
        idx_type col = *col_iter;
        idx_type row = 0;
        for (auto row_iter = rows.cbegin(); row_iter < rows.cend(); row_iter ++)
        {
            auto row_iter_beg = m_row_idx.cbegin() + m_col_idx[col];
            auto row_iter_end = m_row_idx.cbegin() + m_col_idx[col + 1];
            auto pos_iter     = std::equal_range(row_iter_beg, row_iter_end, *row_iter);

            /** 
             * std::equal_range returns a pair of iterators, [first,
             * second] and second = first + 1 if the value to be
             * sought is in the range to be searched.
             */
            if (pos_iter.second > pos_iter.first)
            {
                /** 
                 * Compute the number of non-zero elements in each
                 * column.  Nota bene: I write it to the next column
                 * so I can just compute the column offset using
                 * std::partial_sum
                 */ 

                result.m_col_idx [col_iter - cols.cbegin() + 1] ++;
                result.m_row_idx.push_back (row_iter - rows.cbegin());
                result.m_x.push_back       (*(m_x.cbegin() + std::distance(m_row_idx.cbegin(), pos_iter.first)));
            }
        }
    }

    std::partial_sum(result.m_col_idx.cbegin(), result.m_col_idx.cend(), result.m_col_idx.begin());
    
    /**
     * No need to sort result, the row indices are sorted by
     * construction
     */
    return result;
}
#endif
