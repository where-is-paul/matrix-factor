//-*- mode: c++ -*-
#ifndef HALF_MATRIX_TO_STRING_H
#define HALF_MATRIX_TO_STRING_H

template <class el_type>
std::string half_matrix<el_type>::to_string() const
{
	std::ostringstream os;
	os << "List of Lists Matrix (" << m_n_rows << ", " << m_n_cols << ", " << nnz() << ")" << std::endl;
  
	for (int i = 0; i < n_cols(); i++)
	{
		os << "Column " << i << ":" << std::endl;
		os << "Row      Indices = ["  << m_idx[i] << "]" << std::endl;
		os << "Non-zero Values  = ["  << m_x[i] << "]" << std::endl;
		os << std::endl;
	}

	return os.str();
}

#endif