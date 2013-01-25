// -*- mode: c++ -*-
#ifndef SQUARE_MATRIX_SWAP_STRUCT_H
#define SQUARE_MATRIX_SWAP_STRUCT_H

/*!	\brief A structure containing variables used in pivoting a LIL-C matrix.
	
	Storing these variables in a combined structure reduces memory requirements and bundles together all temporary structures needed during pivoting.
*/
template<class el_type> 
class square_matrix_swap_struct
{
	//---------- useful typedefs (to keep consistent with lilc_matrix) -----------//
	typedef vector<int> idx_vector_type;
	typedef vector<el_type>  elt_vector_type;
	typedef typename idx_vector_type::iterator idx_it;
	typedef typename elt_vector_type::iterator elt_it;
	
public:

	idx_vector_type row_k;	///<Column indices of non-zeros in the new row k.
	idx_vector_type	row_r;	///<Column indices of non-zeros in the new row r.
		
	idx_vector_type col_k_nnzs;	///<Row indices of non-zeros in the new column k.
	elt_vector_type col_k;	///<Non-zero values in the new column k (order dependent on col_k_nnzs).
	idx_vector_type col_r_nnzs;	///<Row indices of non-zeros in the new column r.
	elt_vector_type col_r;	///<Non-zero values in the new column r (order dependent on col_r_nnzs).

	vector<idx_it> swapk;	///<List of indices from row r that will be swapped to row k.
	vector<idx_it> swapr;	///<List of indices from row k that will be swapped to row r.
		
	idx_vector_type all_swaps;	///<Column indices of all swaps done in swapk and swapr.
	

	/*!	\brief Clears all vectors except all_swaps.
	*/
	void clear()
	{
		row_k.clear();
		row_r.clear();

		col_k_nnzs.clear();
		col_k.clear();
		col_r_nnzs.clear();
		col_r.clear();

		swapk.clear();
		swapr.clear();
	}
};

#endif