//-*-mode:c++-*-
#ifndef _SKEW_LILC_MATRIX_LOAD_H_
#define _SKEW_LILC_MATRIX_LOAD_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

template <class el_type>
inline bool skew_readline(std::stringstream& line, int& n_rows, int& n_cols, int& i, int& j, el_type& value)
{
	line >> i >> j >> value;
	i--;
	j--;
	if(i>=0 && j>=0 && i<n_rows&& j< n_cols && i!=j)
		return true;
	else
		return false;
}

template <class el_type>
bool skew_lilc_matrix<el_type>::load(std::string filename)
{
	std::ifstream input(filename.c_str(), std::ios::in);

	if(!input)
		return false;
	
	const int maxBuffersize = 2048;
	char buffer[maxBuffersize];

	bool readsizes = false;

	int n_rows(-1), n_cols(-1), n_nzs(-1), i(-1), j(-1);
	int count = 0;
	el_type value; 

	while (input.getline(buffer, maxBuffersize))
	{
		// skip comments   
		//NOTE An appropriate test should be done on the header to get the skew-symmetry
		if (buffer[0]=='%')
			continue;
		
		std::stringstream line(buffer);
		
		if (!readsizes)
		{
			line >> n_rows >> n_cols >> n_nzs;
			if (n_rows > 0 && n_cols > 0 && n_nzs > 0)
			{
				readsizes = true;
				
				resize(n_rows, n_cols);
				std::fill(first.begin(), first.end(), 0); //a bit of optimization could be used here since resize sets all elem in first to 1
			}
		}
		else
		{ 
			i = -1;
			j = -1;
			if (skew_readline(line, n_rows, n_cols, i, j, value))
			{
				m_idx[j].push_back(i);
				m_x[j].push_back(value);
				count++;
				list[i].push_back(j);
			}
			else
				std::cerr << "Invalid read: " << i << "," << j << "\n";		
		}		
	}
	
	if (count != n_nzs)
		std::cout << "Expected " << n_nzs << " elems but read " << count << "." << std::endl;
	
	nnz_count = count;
	std::cout << "Load succeeded. " << "File " << filename << " was loaded." << std::endl;
	input.close();
	return true;
}

#endif
