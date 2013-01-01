//-*- mode: c++ -*-
#ifndef ULTRIANGULAR_MATRIX_TO_STRING_H
#define ULTRIANGULAR_MATRIX_TO_STRING_H

#include <string>
#include <sstream>

template <class el_type>
std::string ultriangular_matrix<el_type>::to_string() const
{
	std::ostringstream os;
	os << "List of Lists Matrix (" << m_n_rows << ", " << m_n_cols << ", " << nnz() << ")" << std::endl;
  
	for (int i = 0; i < n_cols(); i++)
	{
		int size = m_idx[i].size();
		os << "Column " << i << ":" << std::endl;
		os << "Row      Indices = [" << i;
		if (size > 0)
			os << ", " << m_idx[i];
		os << "]" << std::endl;
		os << "Non-zero Values  = [1";
		if (size > 0)
			os << ", "  << m_x[i];
		os << "]" << std::endl << std::endl;
	}

	return os.str();
}

#endif