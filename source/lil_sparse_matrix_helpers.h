// -*- mode: c++ -*-
#ifndef LIL_SPARSE_MATRIX_HELPERS_H
#define LIL_SPARSE_MATRIX_HELPERS_H

using std::vector;
using std::abs;

typedef vector<int>::iterator idx_it;

const long double eps = 1e-13; //<Machine epsilon

#ifndef DEBUG
#define DEBUG
template<class el_type>
std::ostream& operator<< (std::ostream& os, const vector<el_type>& vec)
{
	if (!vec.empty())
	{
		for (unsigned int index = 0; index < vec.size() - 1; index++)
			os << vec[index] << ", ";

		os << vec[vec.size()-1];
	}
	return os;
}
#endif

/*! \brief Computes the maximum (in absolute value) element of v(curr_nnzs) and it's index.
	\param v the vector whose max element is to be computed.
	\param curr_nnzs a list of indices representing non-zero elements in v.
	\param r the index of the maximum element of v
	\return the absolute value of the max element of v.
*/
template <class el_type>
inline double max(vector<el_type>& v, vector<int>& curr_nnzs, int& r)
{
	double res = 0;
	for (idx_it it = curr_nnzs.begin(), end = curr_nnzs.end(); it != end; it++)
	{
		if (abs(v[*it]) > res)
		{
			res = abs(v[*it]);
			r = *it;
		}
	}
	
	return res;
}

/*! \brief Computes the norm of v(curr_nnzs).
	\param v the vector whose norm is to be computed.
	\param curr_nnzs a list of indices representing non-zero elements in v.
	\param p The norm number.
	\return the norm of v.
*/
template <class el_type>
inline double norm(vector<el_type>& v, vector<int>& curr_nnzs, el_type p = 1)
{
	el_type res = 0;
	for (idx_it it = curr_nnzs.begin(), end = curr_nnzs.end(); it != end; it++)
		res += pow(abs(v[*it]), p);

	return pow(res, 1/p);
}

/*! \brief Functor for comparing elements by value (in decreasing order) instead of by index.
	\param v the vector that contains the values being compared.
*/
template <class el_type>
struct by_value {
	const vector<el_type>& v; 
	by_value(const vector<el_type>& vec) : v(vec) {}
	bool operator()(int const &a, int const &b) const { 
		if (abs(v[a]) == abs(v[b])) return a < b;
		return abs(v[a]) > abs(v[b]);
	}
};

/*! \brief Functor for determining if a variable is below the tolerance given.
    \param v the vector that contains the values being checked.
    \param eps the tolerance given.
*/
template <class el_type>
struct by_tolerance {
  const vector<el_type>& v; 
  double eps;
	by_tolerance(const vector<el_type>& vec, const double& eps) : v(vec), eps(eps) {}
	bool operator()(int const &i) const { 
		return abs(v[i]) < eps;
	}
};

/*! \brief Performs the dual-dropping criteria outlined in Li & Saad (2005).
	\param v the vector that whose elements will be selectively dropped.
	\param curr_nnzs the non-zeros in the vector v.
	\param lfil a parameter to control memory usage. Each column is guarannted to have fewer than lfil elements.
	\param tol a parameter to control agressiveness of dropping. Elements less than tol*norm(v) are dropped.
*/
template <class el_type>
inline void drop_tol(vector<el_type>& v, vector<int>& curr_nnzs, const int& lfil, const double& tol)
{
	//determine dropping tolerance. all elements with value less than tolerance = tol * norm(v) is dropped.
	el_type tolerance = tol*norm(v, curr_nnzs);
	if (tolerance > eps)
	{
		for (idx_it it = curr_nnzs.begin(), end = curr_nnzs.end(); it != end; it++) 
		if (abs(v[*it]) < tolerance) v[*it] = 0;
	}
	
	//sort the remaining elements by value in decreasing order.
	if (lfil < (int) curr_nnzs.size())
	{ //only sort if we can't keep all of them
		by_value<el_type> sorter(v);
		std::sort(curr_nnzs.begin(), curr_nnzs.end(), sorter);
		for (int i = lfil, end = curr_nnzs.size(); i < end ; i++)
			v[curr_nnzs[i]] = 0;
	}
	
	by_tolerance<el_type> is_zero(v, eps);
	curr_nnzs.erase( remove_if(curr_nnzs.begin(), curr_nnzs.end(), is_zero), curr_nnzs.end() );
	curr_nnzs.resize( std::min(lfil, (int) curr_nnzs.size()) );
}

/*! \brief Performs an inplace union of two sorted lists (a and b), removing duplicates in the final list.
	\param a the sorted list to contain the final merged list.
	\param b_start an iterator to the start of b.
	\param b_end an iterator to the end of b.
*/
template <class InputContainer, class InputIterator>
inline void inplace_union(InputContainer& a, InputIterator const& b_start, InputIterator const& b_end)
{
	int mid = a.size(); //store the end of first sorted range

	//copy the second sorted range into the destination vector
	std::copy(b_start, b_end, std::back_inserter(a));

	//perform the in place merge on the two sub-sorted ranges.
	std::inplace_merge(a.begin(), a.begin() + mid, a.end());

	//remove duplicate elements from the sorted vector
	a.erase(std::unique(a.begin(), a.end()), a.end());
}

/*! \brief Performs an inplace union of two unsorted lists (a and b), removing duplicates in the final list.
	\param a the sorted list to contain the final merged list.
	\param b_start an iterator to the start of b.
	\param b_end an iterator to the end of b.
	\param in_set a bitset used to indicate elements present in a and b. Reset to all zeros after union.
*/
template <class InputContainer, class InputIterator>
inline void unordered_inplace_union(InputContainer& a, InputIterator const& b_start, InputIterator const& b_end, vector<bool>& in_set)
{
	for (InputIterator it = a.begin(), end = a.end(); it != end; it++)
		in_set[*it] = true;
	
	for (InputIterator it = b_start; it != b_end; it++)
		if (!in_set[*it])
			a.push_back(*it);
	
	for (InputIterator it = a.begin(), end = a.end(); it != end; it++)
		in_set[*it] = false;
}

inline void safe_swap(vector<int>& curr_nnzs, const int& k, const int& r)
{
	for (idx_it it = curr_nnzs.begin(), end = curr_nnzs.end(); it != end; it++)
	{
		if (*it == k)
			*it = r;
		else if (*it == r)
			*it = k;
	}
}

#endif