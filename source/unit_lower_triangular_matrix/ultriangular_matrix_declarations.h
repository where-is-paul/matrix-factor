// -*- mode: c++ -*-
#ifndef ULTRIANGULAR_MATRIX_DECLARATIONS_H
#define ULTRIANGULAR_MATRIX_DECLARATIONS_H

/*! \brief A list-of-lists (LIL) matrix in column oriented format.

	For convience, the matrix this class represents will be refered to as matrix A.
	In LIL-C format, each column of A (an n*n matrix) is stored as a separate vector. The nonzeros are stored in m_idx while the non-zeros are stored in m_x. Both m_x and m_idx are initialized to a list of n lists. m_idx and m_x are ordered dependent on each other, in that A(m_idx[k][j], k) = m_x[k][j].
	
*/

template <class el_type> 
class ultriangular_matrix : public lil_sparse_matrix<el_type>
{
public:
	
	//-------------- typedefs and inherited variables --------------//
	using lil_sparse_matrix<el_type>::m_idx; // linux compiler gcc need these
	using lil_sparse_matrix<el_type>::m_x;
	using lil_sparse_matrix<el_type>::m_n_rows;
	using lil_sparse_matrix<el_type>::m_n_cols;
	using lil_sparse_matrix<el_type>::n_rows;
	using lil_sparse_matrix<el_type>::n_cols;
	using lil_sparse_matrix<el_type>::nnz;
	using lil_sparse_matrix<el_type>::nnz_count;
	using lil_sparse_matrix<el_type>::list;

	typedef typename lil_sparse_matrix<el_type>::idx_vector_type idx_vector_type; // linux compiler gcc need these
	typedef typename lil_sparse_matrix<el_type>::elt_vector_type elt_vector_type;
	typedef typename lil_sparse_matrix<el_type>::idx_it          idx_it;
	typedef typename lil_sparse_matrix<el_type>::elt_it          elt_it;
	

	idx_vector_type column_first;	///<On iteration k, first[i] gives the number of non-zero elements on col (or row) i of A before A(i, k).
	
	
	/*! \brief Constructor for a column oriented list-of-lists (LIL) matrix. Space for both the values list and the indices list of the matrix is allocated here.
	*/
	ultriangular_matrix (int n_rows = 0, int n_cols = 0): 
	lil_sparse_matrix<el_type> (n_rows, n_cols) 
	{
		m_x.reserve(n_cols);
		m_idx.reserve(n_cols);
	}
	
	//----Matrix referencing/filling----//
	
	/*! \brief Resizes the matrix. For use in preallocating space before factorization begins.
		\param n_rows the number of rows in the resized matrix.
		\param n_cols the number of cols in the resized matrix.
	*/
	void resize(int n_rows, int n_cols)
	{
		lil_sparse_matrix<el_type>::resize(n_rows, n_cols); // call the function resize in base class

		column_first.resize(n_cols, 0);
	}


	//------Helpers------//
	/*! \brief Ensures the two invariants observed by L.column_first and L.m_idx held.
		
		\invariant
		If this matrix is a lower triangular factor of another matrix:
			-# On iteration k, first[i] will give the number of non-zero elements on col i of A before A(k, i).
			-# On iteration k, list[i][ first[i] ] will contain the first element below or on index k of column i of A.
		
		\invariant
		If this matrix is the matrix to be factored:
			-# On iteration k, first[i] will give the number of non-zero elements on row i of A before A(i, k).
			-# On iteration k, list[i][ first[i] ] will contain the first element right of or on index k of row i of A.
			
		\param j the column of con.
		\param k the iteration number.
		\param con the container to be swapped.
		\param update_list boolean indicating whether list or m_x/m_idx should be updated.
	*/
	inline void ensure_invariant(const int& j, const int& k, idx_vector_type& con)
	{
		int offset = column_first[j];
		if (con[offset] == k)
			return;
		
		int i, min = offset;
		for (i = offset; i < (int) con.size(); i++)
		{
			if (con[i] == k)
			{
				min = i;
				break;
			}
		}
		
		std::swap(con[offset], con[min]);
		std::swap(m_x[j][offset], m_x[j][min]);
	}
	
	/*! \brief Updates L.first for iteration k.
		\param k current iteration index.	
	*/
	inline void advance_column(const int& k)
	{
		for (auto it = list[k].begin(); it != list[k].end(); it++)
		{	
			ensure_invariant(*it, k, m_idx[*it]); //make sure next element is good before we increment.
			column_first[*it]++; //should have ensured invariant now
		}
	}
	
	//----IO Functions----//
	
	/*! \brief Returns a string representation of A, with each column and its corresponding indices & non-zero values printed.
		\return A string representation of this matrix.
	*/
	std::string to_string () const;
	
	/*! \brief Saves a matrix in matrix market format.
		\param filename the filename of the matrix to be saved. All matrices saved are in matrix market format (.mtx).
		\param sym flags whether the matrix is symmetric or not.
	*/
	bool save(std::string filename);

};

//------------------ include files for class functions -------------------//

#include "ultriangular_matrix_save.h"
#include "ultriangular_matrix_to_string.h"

#endif