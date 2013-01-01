//-*-mode:c++-*-
#ifndef ULTRIANGULAR_MATRIX_SAVE_H
#define ULTRIANGULAR_MATRIX_SAVE_H

inline void ultriangular_put_header(std::string& header)
{
	header= "%%MatrixMarket matrix coordinate real ";
	header += "general";
}

template <class el_type>
bool ultriangular_matrix<el_type>::save(std::string filename)
{
	std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
	if(!out)
	return false;

	out.flags(std::ios_base::scientific);
	out.precision(16);
	std::string header; 
	ultriangular_put_header(header); 

	out << header << std::endl; 
	out << n_rows() << " " << n_cols() << " " << nnz() << "\n";

	for(int i = 0; i < n_cols(); i++) {
		out << i+1 << " " << i+1 << " " << 1 << "\n";
		for(unsigned int j = 0; j < m_idx[i].size(); j++) {
			out << m_idx[i][j]+1 << " " << i+1 << " " << m_x[i][j] << "\n";
		}
	}
	
	out.close();
	return true;
}

#endif