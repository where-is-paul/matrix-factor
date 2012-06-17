//-*-mode:c++-*-

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: compress_sum (const triplet_matrix<idx_type, el_type>& triplet)
{
    const std::vector<idx_type>& t_row_idx = triplet.row_idx();
    const std::vector<idx_type>& t_col_idx = triplet.col_idx();
    const std::vector<el_type>& t_data    = triplet.nz_vals();

    m_n_rows = triplet.n_rows();
    m_n_cols = triplet.n_cols();
  
    // Compute the number of non-zero entries in each column
    std::vector<idx_type> work(m_n_cols);
    std::for_each (t_col_idx.begin(), t_col_idx.end(), [&work](idx_type col) {work[col] ++;});
  
    // Compute column indices as the cumsum of the non-zero entries in each column
    m_col_idx.resize (n_cols() + 1);
    m_col_idx[0] = 0;
    std::partial_sum (work.begin(), work.end(), m_col_idx.begin() + 1);
    // Overwrite work to contain the start indices of the rows in each column  
    std::copy (m_col_idx.begin(), m_col_idx.end() - 1, work.begin());

    idx_type nnz = triplet.nnz();
    m_row_idx.resize (nnz);
    m_x.resize (nnz);

    for (idx_type elem = 0; elem < nnz; elem ++)
    {
        idx_type& index = work[t_col_idx[elem]];
        m_row_idx[index] = t_row_idx [elem];
        m_x [index] = t_data    [elem];
        index ++;
    }

}
