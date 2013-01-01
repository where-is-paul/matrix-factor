//-*-mode:c++-*-
#ifndef HALF_MATRIX_SAVE_H
#define HALF_MATRIX_SAVE_H

template <class el_type>
bool half_matrix<el_type>::save(std::string filename)
{
	std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
	if(!out)
		return false;

	out.flags(std::ios_base::scientific);
	out.precision(16);
	std::string header;
	put_header(header);

	out << header << std::endl; 

	std::stringstream fout;

	int non_zeros = 0;
	for(int i = 0; i < n_cols(); i++) {
		for(unsigned int j = 0; j < m_idx[i].size(); j++) {
			fout << m_idx[i][j]+1 << " " << i+1 << " " << m_x[i][j] << "\n";
			non_zeros++;
		}
	}
	
	out << n_rows() << " " << n_cols() << " " << non_zeros << "\n";
	out << fout.str();
	out.close();
	return true;
}

#endif