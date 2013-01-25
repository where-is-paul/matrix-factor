// -*- mode: c++ -*-
#ifndef SQUARE_MATRIX_DECLARATIONS_H
#define SQUARE_MATRIX_DECLARATIONS_H

#include <iterator>
#include "square_matrix_swap_struct.h"
#include "../unit_lower_triangular_matrix/ultriangular_matrix.h"

/*! \brief A list-of-lists (LIL) matrix in column oriented format.

	For convience, the matrix this class represents will be refered to as matrix A.
	In LIL-C format, each column of A (an n*n matrix) is stored as a separate vector. The nonzeros are stored in m_idx while the non-zeros are stored in m_x. Both m_x and m_idx are initialized to a list of n lists. m_idx and m_x are ordered dependent on each other, in that A(m_idx[k][j], k) = m_x[k][j].
	
*/
template <class el_type> 
class square_matrix : public lil_sparse_matrix<el_type>
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


	idx_vector_type list_first;	///<On iteration k, first[i] gives the number of non-zero elements on col (or row) i of A before A(i, k).

	elt_vector_type S; ///<A diagonal scaling matrix S such that SAS will be equilibriated in the max-norm (i.e. every row/column has norm 1). S is constructed after running the sym_equil() function, after which SAS will be stored in place of A.
	int sign;	/// 1 for symmetric, -1 for skew-symmetric


	/*! \brief Constructor for a column oriented list-of-lists (LIL) matrix. Space for both the values list and the indices list of the matrix is allocated here.
	*/
	square_matrix(int n_rows = 0, int n_cols = 0) : lil_sparse_matrix<el_type> (n_rows, n_cols)
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
	inline void coeffRef(const int& i, const int& j, idx_it& idx_iterator, elt_it& elt_iterator)
	{
		for (unsigned int k = 0; k < m_idx[j].size(); k++)
		{
			if (m_idx[j][k] == i)
			{
				idx_iterator = m_idx[j].begin() + k;
				elt_iterator = m_x[j].begin() + k;
				break;
			}
		}
	}

	/*! \brief Resizes the matrix. For use in preallocating space before factorization begins.
		\param n_rows the number of rows in the resized matrix.
		\param n_cols the number of cols in the resized matrix.
	*/
	void resize(int n_rows, int n_cols)
	{
		lil_sparse_matrix<el_type>::resize(n_rows, n_cols); // call the function resize in base class

		list_first.resize(n_cols, 0);
		S.resize(n_cols, 0);
	}
	
	//-----Reorderings/Rescalings------//

	/*!	\brief The symmetric matrix A is equilibrated and the symmetric equilibrated matrix SAS is stored in A, where S is a diagonal scaling matrix. 
		This algorithm is based on the one outlined in "Equilibration of Symmetric Matrices in the Max-Norm" by Bunch (1971).
	*/
	void equilibrate();

	/*!	\brief Returns a pseudo-peripheral root of A. This is essentially many chained breadth-first searchs across the graph of A (where A is viewed as an adjacency matrix).

		\param s contains the initial node to seed the algorithm. A pseudo-peripheral root of A is stored in s at the end of the algorithm.
	*/
	inline void find_root(int& s, vector<idx_vector_type>& adjacency_list, vector<int>& deg);

	/*!	\brief Returns the next level set given the current level set of A. This is essentially all neighbours of the currently enqueued nodes in breath-first search.
		
		\param lvl_set the current level set (a list of nodes).
		\param visited all previously visited nodes.
	*/
	//inline bool find_level_set(vector<int>& lvl_set, vector<bool>& visited);

	/*!	\brief Returns a Reverse Cuthill-McKee ordering of the matrix A (stored in perm). 
		
		The implementation is based on the general algorithm outlined in A detailed description of this function as well as all its subfunctions can be found in "Computer Solution of Large Sparse Positive Definite Systems" by George and Liu (1981).
		\param perm An empty permutation vector (filled on function completion).
	*/
	inline void rcm(vector<int>& perm);

	/*! \brief Given a permutation vector perm, A is permuted to P'AP, where P is the permutation matrix associated with perm. 
		\param perm the permutation vector.
	*/
	inline void permute(vector<int>& perm);

	//------Helpers------//
	/*! \brief Ensures the two invariants observed by A.list_first and A.list are held.
		
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
		int offset = list_first[j];
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
	}

	/*! \brief Performs a symmetric permutation between row/col k & r of A.
		\param s a struct containing temporary variables needed during pivoting.
		\param in_set a bitset needed for unordered unions during pivoting.
		\param L the lower triangular factor of A.
		\param k index of row/col k.
		\param r index of row/col r.
	*/
	inline void pivot(square_matrix_swap_struct<el_type> s, ultriangular_matrix<el_type>& L, const int& k, const int& r);
	
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
	bool save(std::string filename1, std::string filename2);

		/*!	\brief Finds the degree of a node.
		\param index the index of the node.
		\return the degree of the nodex indexed by index.
	*/
	virtual int degree(const int& index) const = 0;

	virtual void put_header(std::string& header) = 0;

	virtual bool readline (char*& line, int& n_rows, int& n_cols, int& i, int& j, el_type& value) = 0;

};

//------------------ include files for class functions -------------------//

#include "square_matrix_equilibrate.h"
//#include "square_matrix_find_level_set.h"
#include "square_matrix_find_root.h"
#include "square_matrix_rcm.h"
#include "square_matrix_load.h"
#include "square_matrix_save.h"
#include "square_matrix_to_string.h"
#include "square_matrix_permute.h"
#include "square_matrix_pivot.h"

#endif