//-*-mode:c++-*-
#ifndef _SKEW_LILC_MATRIX_SAVE_H_
#define _SKEW_LILC_MATRIX_SAVE_H_

inline void skew_put_header(std::string& header, bool sym = false)
{
	header= "%%MatrixMarket matrix coordinate real ";
	if (sym)
		header += "skew-symmetric"; //maybe change later to have skew-symmetric/complex/blah as options
	else
		header += "general";
}

template <class el_type>
bool skew_lilc_matrix<el_type>::save(std::string filename, bool sym = false)
{
	std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
	if(!out)
		return false;

	out.flags(std::ios_base::scientific);
	out.precision(16);
	std::string header; 
	skew_put_header(header, sym); 

	out << header << std::endl; 
	out << n_rows() << " " << n_cols() << " " << nnz() << "\n";

	for(int j = 0; j < n_cols(); j++)
		for(unsigned int i = 0; i < m_idx[j].size(); i++)
			out << m_idx[j][i]+1 << " " << j+1 << " " << m_x[j][i] << "\n";
	
	out.close();
	return true;
}

#endif