//-*-mode:c++-*-
#ifndef HALF_MATRIX_SAVE_H
#define HALF_MATRIX_SAVE_H

template <class el_type>
bool square_matrix<el_type>::save(std::string filename1, std::string filename2)
{
	int ncols = n_cols();

	//------------------------start save A------------------------------
	std::ofstream out(filename1.c_str(), std::ios::out | std::ios::binary);
	if (!out)
		return false;
	std::stringstream fout;
	std::string header;

	//out.flags(std::ios_base::scientific);
	out.precision(16);
	//fout.flags(std::ios_base::scientific);
	fout.precision(16);
	
	put_header(header);
	out << header << std::endl;

	int non_zeros = 0;
	for (int i = 0; i < ncols; i++)
	{
		for (int j = 0, size = m_idx[i].size(); j < size; j++)
		{
			fout << m_idx[i][j]+1 << " " << i+1 << " " << m_x[i][j] << "\n";
			non_zeros++;
		}
	}
	
	out << ncols << " " << ncols << " " << non_zeros << "\n";
	out << fout.str();
	out.close();
	//--------------------------end save A------------------------------

	//------------------------start save S------------------------------
	out.open(filename2.c_str(), std::ios::out | std::ios::binary);
	if (!out)
		return false;

	//out.flags(std::ios_base::scientific);
	out.precision(16);
	
	header = "%%MatrixMarket matrix coordinate real general";
	out << header << std::endl;
	out << ncols << " " << ncols << " " << ncols << "\n";

	for (int i = 0; i < ncols; i++)
		out << i+1 << " " << i+1 << " " << S[i] << "\n";
	
	out.close();
	//--------------------------end save S------------------------------

	return true;
}

#endif