//-*-mode:c++-*-
#ifndef _SKEW_BLOCK_DIAG_MATRIX_SAVE_H_
#define _SKEW_BLOCK_DIAG_MATRIX_SAVE_H_

template <class el_type>
bool skew_block_diag_matrix<el_type>::save(std::string filename) const
{
	std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
	if (!out)
		return false;

	//out.flags(std::ios_base::scientific);
	out.precision(16);
	std::string header;
	
	header = "%%MatrixMarket matrix coordinate ";
	header += "real skew-symmetric"; //maybe change later to have general/complex/blah as options

	out << header << std::endl; 
	out << n_rows() << " " << n_cols() << " " << nnz() << "\n";

	for(int i = 1; i < n_cols(); i=i+2)
		out << i+1 << " " << i << " " << subdiag[i/2] << "\n";
	
	out.close();
	return true;
}

#endif
