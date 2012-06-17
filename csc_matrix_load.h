//-*-mode:c++-*-
#ifndef _CSC_MATRIX_LOAD_H_
#define _CSC_MATRIX_LOAD_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

template <class idx_type, class el_type>
inline bool readline (std::stringstream& line, idx_type& n_rows, idx_type& n_cols, idx_type& i, idx_type& j, el_type& value)
{
    line >> i >> j >> value;
    i--;
    j--;
    if(i>=0 && j>=0 && i<n_rows&& j< n_cols)
    {
      return true; 
    }
    else
      return false;
}
  
template <class idx_type, class el_type>
bool csc_matrix<idx_type, el_type> :: load (std::string filename)
{
  std::ifstream input(filename.c_str(), std::ios::in);

  if(!input) return false;
  
  const int maxBuffersize = 2048;
  char buffer[maxBuffersize];
  
  bool readsizes = false;
  
  idx_type n_rows(-1), n_cols(-1), n_nzs(-1), prev_j(0);
  int count = 0, prev_count = 0; 
  while(input.getline(buffer, maxBuffersize))
  {
    // skip comments   
    //NOTE An appropriate test should be done on the header to get the symmetry
    if(buffer[0]=='%')
      continue;
    
    std::stringstream line(buffer);
    
    if(!readsizes)
    {
      line >> n_rows >> n_cols >> n_nzs;
      if(n_rows > 0 && n_cols > 0 && n_nzs > 0) 
      {
        readsizes = true;
        //std::cout << "Sizes: " << n_rows << ", " << n_cols << ", " << n_nzs << "\n";
		resize(n_rows, n_cols, n_nzs);
      }
    }
    else
    { 
      idx_type i(-1), j(-1);
      el_type value; 
      if( readline(line, n_rows, n_cols, i, j, value) ) 
      {
		if (prev_j != j)
		{	
			//std::cout << prev_j << " transitioning to " << j << " on element number: " << count << std::endl;
			for (idx_type k = prev_j; k < j; k++) m_col_idx[k] = prev_count;
			m_col_idx[j] = count;
			prev_count = count;
			prev_j = j;
		}
		
		m_row_idx[count] = i;
        m_x[count] = value;
		++count;
      }
      else 
        std::cerr << "Invalid read: " << i << "," << j << "\n";		
    }
	
	m_col_idx[n_cols] = n_nzs;
  }

  if(count!=n_nzs)
    std::cerr << count << " elements read but expected " << n_nzs << "elements. \n";
  
  std::cout << "Load succeeded. " << "File " << filename << " was loaded." << std::endl;
  input.close();
  return true;
}


#endif
