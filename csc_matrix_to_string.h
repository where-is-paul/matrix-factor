//-*- mode: c++ -*-
#ifndef _CSC_MATRIX_TO_STRING_H_
#define _CSC_MATRIX_TO_STRING_H_

#include <string>
#include <sstream>

template <class idx_type, class el_type>
std::string csc_matrix<idx_type, el_type> :: to_string() const
{
  std::ostringstream os;
  os << "Compressed Sparse Column Matrix (" << m_n_rows << ", " << m_n_cols << ", " << nnz() << ")" << std::endl;
  
  os << "Column   Indices = "  << m_col_idx << std::endl;
  os << "Row      Indices = "  << m_row_idx << std::endl;
  os << "Non-zero Values  = "  << m_x       << std::endl;
  return os.str();
}

#endif // _CSC_MATRIX_TO_STRING_H_
