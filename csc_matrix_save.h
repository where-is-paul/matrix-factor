//-*-mode:c++-*-
#ifndef _CSC_MATRIX_SAVE_H_
#define _CSC_MATRIX_SAVE_H_

inline void put_header(std::string& header)
{
	header= "%%MatrixMarket matrix coordinate ";
	header += "real general"; //maybe change later to have symmetric/complex/blah as options
}
  
template <class idx_type, class el_type>
bool csc_matrix<idx_type, el_type> :: save(std::string filename)
{
  std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
  if(!out)
    return false;

  out.flags(std::ios_base::scientific);
  out.precision(10);
  std::string header; 
  put_header(header); 
  
  out << header << std::endl; 
  out << n_rows() << " " << n_cols() << " " << nnz() << "\n";
  
  for(idx_type i = 0; i < (idx_type) m_col_idx.size()-1; i++)
    for(idx_type j = m_col_idx[i]; j < (idx_type) m_col_idx[i+1]; j++) {
		out << i+1 << " " << m_row_idx[j]+1 << " " << m_x[j] << "\n";
    }
	
  out.close();
  return true;
}

#endif