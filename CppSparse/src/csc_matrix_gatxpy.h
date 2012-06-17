//-*-mode:c++-*-
#ifndef _CSC_MATRIX_GATXPY_H_
#define _CSC_MATRIX_GATXPY_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: gatxpy (const std::vector<el_type>& x, std::vector<el_type>& y) const
{
  for (idx_type col = 0; col < m_n_cols; col ++)
  {
      for (idx_type offset = m_col_idx [col]; offset < m_col_idx [col + 1]; offset ++)
      {
          y [col] += m_x [offset] * x[m_row_idx[offset]];
      }
  }
}

#endif // _CSC_MATRIX_GATXPY_H_
