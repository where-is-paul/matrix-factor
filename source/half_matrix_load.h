//-*-mode:c++-*-
#ifndef HALF_MATRIX_LOAD_H
#define HALF_MATRIX_LOAD_H

template <class el_type>
bool half_matrix<el_type>::load(std::string filename)
{	
	std::ifstream input(filename.c_str(), std::ios::in);
	std::stringstream line;
		
	if(!input)
		return false;
	
	const int maxBuffersize = 2048;
	char buffer[maxBuffersize];

	bool readsizes = false;

	int n_rows(-1), n_cols(-1), n_nzs(-1), i(-1), j(-1);
	int readcount = 0, nnzcount = 0;
	el_type value; 

	while(input.getline(buffer, maxBuffersize))
	{
		// skip comments   
		//NOTE An appropriate test should be done on the header to get the symmetry
		if (buffer[0]=='%')
			continue;
		
		line << buffer;
		
		if (!readsizes)
		{
			line >> n_rows >> n_cols >> n_nzs;
			if (n_rows > 0 && n_cols > 0 && n_nzs > 0) 
			{
				readsizes = true;
				resize(n_rows, n_cols);
			}
		}
		else
		{ 
			i = -1;
			j = -1;
			if (readline(line, n_rows, n_cols, i, j, value)) 
			{
				m_idx[j].push_back(i);
				m_x[j].push_back(value);
				readcount++;

				if (i != j)
				{
					list[i].push_back(j);
					nnzcount = nnzcount + 2;
				}
				else
					nnzcount++;
			}
			else 
				std::cerr << "Invalid read: " << i << "," << j << "\n";		
		}
		
		line.clear();
	}
	
	if (readcount != n_nzs)
		std::cout << "Expected " << n_nzs << " elems but read " << readcount << "." << std::endl;
	
	nnz_count = nnzcount;
	std::cout << "Load succeeded. " << "File " << filename << " was loaded." << std::endl;
	input.close();
	return true;
}

#endif