// -*- mode: c++ -*-
#ifndef HALF_MATRIX_DECLARATIONS_H
#define HALF_MATRIX_DECLARATIONS_H

#include <algorithm>
#include <cassert>
#include <iterator>
#include "symmetry_swap_struct.h"

/*! \brief A list-of-lists (LIL) matrix in column oriented format.

	For convience, the matrix this class represents will be refered to as matrix A.
	In LIL-C format, each column of A (an n*n matrix) is stored as a separate vector. The nonzeros are stored in m_idx while the non-zeros are stored in m_x. Both m_x and m_idx are initialized to a list of n lists. m_idx and m_x are ordered dependent on each other, in that A(m_idx[k][j], k) = m_x[k][j].
	
*/
template <class el_type> 
class half_matrix : public lil_sparse_matrix<el_type>
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

	// using lil_sparse_matrix<el_type>::idx_vector_type;
	// using lil_sparse_matrix<el_type>::elt_vector_type;
	// using lil_sparse_matrix<el_type>::idx_it;
	// using lil_sparse_matrix<el_type>::elt_it;

	typedef typename lil_sparse_matrix<el_type>::idx_vector_type idx_vector_type; // linux compiler gcc need these
	typedef typename lil_sparse_matrix<el_type>::elt_vector_type elt_vector_type;
	typedef typename lil_sparse_matrix<el_type>::idx_it idx_it;
	typedef typename lil_sparse_matrix<el_type>::elt_it elt_it;

	// typedef typename idx_vector_type::iterator idx_it;
	// typedef typename elt_vector_type::iterator elt_it;

	vector<idx_vector_type> list;	///<A list of linked lists that gives the non-zero elements in each row of A. Since at any time we may swap between two rows, we require linked lists for each row of A.
	idx_vector_type list_first;	///<On iteration k, first[i] gives the number of non-zero elements on col (or row) i of A before A(i, k).

	block_diag_matrix<el_type> S; ///<A diagonal scaling matrix S such that SAS will be equilibriated in the max-norm (i.e. every row/column has norm 1). S is constructed after running the sym_equil() function, after which SAS will be stored in place of A.
	
public:
	
	/*! \brief Constructor for a column oriented list-of-lists (LIL) matrix. Space for both the values list and the indices list of the matrix is allocated here.
	*/
	half_matrix(int n_rows = 0, int n_cols = 0) : lil_sparse_matrix<el_type> (n_rows, n_cols)
	{
		m_x.reserve(n_cols);
		m_idx.reserve(n_cols);
	}

	/*! \brief Finds the index/value pointers to (i,j)th coefficient of the matrix.
		\param i the row of the (i,j)th element (zero-indexed).
		\param j the col of the (i,j)th element (zero-indexed).
		\param its a pair of pointers, one for the index of the found element, and the other for the value of the element. If the element is not found, the pointers point to the end of column j.
		
		\return True if (i,j)th element is nonzero, false otherwise. 
	*/
	inline bool coeffRef(const int& i, const int& j, std::pair<idx_it, elt_it>& its)
	{	
		for (unsigned int k = 0; k < m_idx[j].size(); k++)
		{
			if (m_idx[j][k] == i)
			{
				its = make_pair(m_idx[j].begin() + k, m_x[j].begin() + k);
				return true;
			}
		}
		
		its = make_pair(m_idx[j].end(), m_x[j].end());
		return false;
	}

	/*! \brief Resizes the matrix. For use in preallocating space before factorization begins.
		\param n_rows the number of rows in the resized matrix.
		\param n_cols the number of cols in the resized matrix.
	*/
	void resize(int n_rows, int n_cols)
	{
		m_n_rows = n_rows;
		m_n_cols = n_cols;
		m_x.resize(n_cols);
		m_idx.resize(n_cols);
		
		list_first.resize(n_cols, 0);
		list.resize(n_cols);
		
		S.resize(n_cols);
	}
	
	//-----Reorderings/Rescalings------//



	/*!	\brief Returns a pseudo-peripheral root of A. This is essentially many chained breadth-first searchs across the graph of A (where A is viewed as an adjacency matrix).

		\param s contains the initial node to seed the algorithm. A pseudo-peripheral root of A is stored in s at the end of the algorithm.
	*/
	inline void find_root(int& s);

	/*!	\brief Returns the next level set given the current level set of A. This is essentially all neighbours of the currently enqueued nodes in breath-first search.
		
		\param lvl_set the current level set (a list of nodes).
		\param visited all previously visited nodes.
	*/
	inline bool find_level_set(vector<int>& lvl_set, vector<bool>& visited);

	/*!	\brief Returns a Reverse Cuthill-McKee ordering of the matrix A (stored in perm). 
		
		The implementation is based on the general algorithm outlined in A detailed description of this function as well as all its subfunctions can be found in "Computer Solution of Large Sparse Positive Definite Systems" by George and Liu (1981).
		\param perm An empty permutation vector (filled on function completion).
	*/
	void rcm(vector<int>& perm);


	//----IO Functions----//
	
	/*! \brief Returns a string representation of A, with each column and its corresponding indices & non-zero values printed.
		\return A string representation of this matrix.
	*/
	std::string to_string () const;

	/*! \brief Loads a matrix in matrix market format.
		\param filename the filename of the matrix to be loaded. Must be in matrix market format (.mtx).
	*/
	bool load(std::string filename);

	/*! \brief Saves a matrix in matrix market format.
		\param filename the filename of the matrix to be saved. All matrices saved are in matrix market format (.mtx).
	*/
	bool save(std::string filename);

		/*!	\brief Finds the degree of a node.
		\param index the index of the node.
		\return the degree of the nodex indexed by index.
	*/
	virtual int degree(const int& index) const = 0;

	virtual void put_header(std::string& header) = 0;

	virtual bool readline (std::stringstream& line, int& n_rows, int& n_cols, int& i, int& j, el_type& value) = 0;

};

//------------------ include files for class functions -------------------//

#include "half_matrix_find_level_set.h"
#include "half_matrix_find_root.h"
#include "half_matrix_matrix_rcm.h"
#include "half_matrix_load.h"
#include "half_matrix_save.h"
#include "half_matrix_to_string.h"

#endif