//-*- mode: c++ -*-
#ifndef _CSC_MATRIX_SORT_H_
#define _CSC_MATRIX_SORT_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: sort ()
{
    // This is just the transpose step
    //  - Compute the number of non-zero columns in each row
    //  - Compute the cumulative sum to create the col_idx vector
    //  - Loop over the entries of the original matrix to create the transpose

    // Work needs to be used twice: to populate the non-zero indices
    // in A' followed by the non-zero indices in A
    std::vector<idx_type> work (std::max(n_rows(), n_cols()));

    std::for_each(m_row_idx.cbegin(), m_row_idx.cend(), [&work](idx_type row) {work[row] ++;});
    csc_matrix<idx_type, el_type> trans(n_cols(), n_rows(), nnz());
    trans.m_row_idx.resize(nnz());
    trans.m_x.resize(nnz());
    trans.m_col_idx[0] = 0;
    std::partial_sum(work.cbegin(), work.cbegin() + n_rows(), trans.m_col_idx.begin() + 1);
    std::copy(trans.m_col_idx.cbegin(), trans.m_col_idx.cbegin() + n_rows(), work.begin());
    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {
            idx_type& index        = work[m_row_idx[offset]];
            trans.m_row_idx[index] = col;
            trans.m_x      [index] = m_x[offset];
            index ++;
        }
    }

    // As said in the text book, we already know the number of
    // non-zero rows in each column of A, so we can just populate the
    // entries without counting them

    std::copy (m_col_idx.begin(), m_col_idx.begin() + n_cols(), work.begin());

    for (idx_type row = 0; row < m_n_rows; row ++)
    {
        for (idx_type offset = trans.m_col_idx[row]; offset < trans.m_col_idx[row + 1]; offset ++)
        {
            idx_type& index   = work[trans.m_row_idx[offset]];
            m_row_idx [index] = row;
            m_x       [index] = trans.m_x[offset];
            index ++;
        }
    }
    
}
#endif
