//-*-mode:c++-*-
#ifndef _TRIPLET_MATRIX_LOAD_H_
#define _TRIPLET_MATRIX_LOAD_H_

#include <fstream>
#include <string>

template <class idx_type, class el_type>
void triplet_matrix<idx_type, el_type> :: load (std::istream& stream)
{
  idx_type row, col;
  el_type  nzval;

  while (stream >> row >> col >> nzval)
  {
    push_back (row, col, nzval);
  }
}


template<class idx_type, class el_type>
triplet_matrix<idx_type, el_type> :: triplet_matrix (std::istream& stream) : abstract_sparse_matrix<idx_type, el_type> (0, 0, 10)
{
  load (stream);
}

template<class idx_type, class el_type>
triplet_matrix<idx_type, el_type> :: triplet_matrix (const std::string& filename) : abstract_sparse_matrix<idx_type, el_type> (0, 0, 10)
{
  std::ifstream file (filename.c_str());
  load (file);
  file.close();
}

#endif
