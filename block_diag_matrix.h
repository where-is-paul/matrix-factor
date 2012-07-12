// -*- mode: c++ -*-
#ifndef _BLOCK_DIAG_MATRIX_H_
#define _BLOCK_DIAG_MATRIX_H_

#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>

#ifndef VECTOR_SHIFT
#define VECTOR_SHIFT
template<class el_type>
std::ostream& operator<< (std::ostream& os, const std::vector<el_type>& vec)
{
	os << "[";
	if (!vec.empty())
	{
		for (typename std::vector<el_type>::size_type index = 0; index < vec.size() - 1; index ++)
		{
			os << vec[index] << ", ";
		}

		os << vec[vec.size()-1];
	}
	os << "]";
	return os;
}
#endif


/*! \brief A quick implementation of a diagonal matrix with 1x1 and 2x2 blocks. */
template<class el_type>
class block_diag_matrix
{
public:

	typedef std::unordered_map<int, el_type> int_elt_map;
	typedef std::vector<el_type>  elt_vector_type;
	
	/*! Allows outputting the contents of the matrix via << operators. */
	friend std::ostream & operator<<(std::ostream& os, const block_diag_matrix& D) 
	{
		os << D.to_string();
		return os;
	};
	
	int m_n_size, nnz_count;
	elt_vector_type main_diag;
	int_elt_map off_diag;
	
	block_diag_matrix (int n_rows = 0, int n_cols = 0) : m_n_size(n_rows)
	{
		nnz_count = n_rows;
		main_diag.resize(n_rows);
	}
	
	void resize(int n) {
		m_n_size = n;
		main_diag.resize(n);
		nnz_count = n;
	}
	
	/*! \return Number of rows in the matrix. */
	const int n_rows() const
	{
		return m_n_size;
	}

	/*! \return Number of cols in the matrix. */
	const int n_cols() const
	{
		return m_n_size;
	}

	/*! \return Number of nonzeros in the matrix. */
	const int nnz() const 
	{
		return nnz_count;
	};
	
	el_type& operator[](int i) {
		return main_diag.at(i);
	}

	el_type& off_diagonal(int i) {
		if (!off_diag.count(i)) {
			off_diag.insert(std::make_pair(i, 0));
			nnz_count++;
		}
		
		return off_diag.at(i);
	}
	
	const int block_size(int i) const {
		if (off_diag.count(i)) {
			return 2;
		} else if (off_diag.count(i-1)) {
			return -2;
		} else {
			return 1;
		}
	}
	
	/*! \return A string reprepsentation of this matrix.
	*/
	std::string to_string () const;
	
	/*! \param filename the filename of the matrix to be saved. All matrices saved are in matrix market format (.mtx).
	*/
	bool save(std::string filename);
	
	~block_diag_matrix()
	{
	}
};

#include "block_diag_matrix_to_string.h"
#include "block_diag_matrix_save.h"

#endif 