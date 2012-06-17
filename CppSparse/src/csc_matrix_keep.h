//-*- mode: c++ -*-
#ifndef _CSC_MATRIX_KEEP_H_
#define _CSC_MATRIX_KEEP_H_

template<class idx_type, class el_type>
template<class predicate>
void csc_matrix<idx_type, el_type> :: keep (const predicate& pred)
{
  idx_type nnz = 0;

  for (idx_type col = 0; col < m_n_cols; col ++)
  {
      idx_type col_start = nnz;
      for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col+1]; offset ++)
      {
          idx_type row = m_row_idx [offset];
          el_type  val = m_x [offset];
          if (pred (row, col, val))
          {
              m_row_idx [nnz] = row;
              m_x       [nnz] = val;
              nnz ++;
          }
      }
      m_col_idx[col] = col_start;
  }
  
  m_col_idx.back() = nnz;
  m_row_idx.resize(nnz);
  m_x.resize(nnz);
}
#endif // _CSC_MATRIX_KEEP_H_
