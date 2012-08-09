template<class el_type> 
class swap_struct
{
	typedef vector<int> idx_vector_type;
	typedef vector<el_type>  elt_vector_type;
	typedef typename idx_vector_type::iterator idx_it;
	typedef typename elt_vector_type::iterator elt_it;
	
	public:
		vector<idx_it> swapk;
		vector<idx_it> swapr;
		
		idx_vector_type all_swaps;
		idx_vector_type col_k_nnzs;
		idx_vector_type col_r_nnzs;
		
		elt_vector_type col_k;
		elt_vector_type col_r;

		void swap_clear() {
			swapk.clear();
			swapr.clear();
			all_swaps.clear();
		}
		
		void col_clear() {
			col_k.clear();
			col_k_nnzs.clear();
			col_r.clear();
			col_r_nnzs.clear();
		}
};