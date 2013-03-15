//-*- mode: c++ -*-
#ifndef SQUARE_MATRIX_EQUILIBRATE_H
#define SQUARE_MATRIX_EQUILIBRATE_H

using std::abs;

template<class el_type>
void square_matrix<el_type>::equilibrate()
{
	//find termination points for loops with binary search later.
	int i, ncols = n_cols();
	idx_it idx_iterator;
	elt_it elt_iterator;
	vector<elt_it> temp(ncols);

	for (i = 0; i < ncols; i++)
	{
		//assumes diag elem is always in 0th pos. if possible.
		if (!m_idx[i].empty() && m_idx[i][0] == i)
			S[i] = sqrt(abs(m_x[i][0]));
		
		//assumes indices are ordered. since this procedure is run
		//before factorization pivots matrix, this is a fair assumption
		//for most matrix market matrices.
		for (auto it = list[i].begin(); it != list[i].end(); it++)
		{
			//if (S[*it] < eps) continue;
			coeffRef(i, *it, idx_iterator, elt_iterator);
			temp[*it] = elt_iterator;
			S[i] = std::max(S[i], abs(*elt_iterator));
		}
		
		//S[i] > 0 since its the square root of a +ve number
		if (S[i] > eps)
		{ 
			for (auto it = list[i].begin(); it != list[i].end(); it++)
				//can use bin. search on coeff since no reordering is done yet.
				*temp[*it] /= S[i];

			if (!m_idx[i].empty() && (m_idx[i][0] == i))
				m_x[i][0] /= S[i];
			for (auto it = m_x[i].begin(); it != m_x[i].end(); it++)
				*it /= S[i];
		} 
	}
	
	for (i = 0; i < ncols; i++)
	{
		if (S[i] < eps)
		{
			for (auto it = m_x[i].begin(); it != m_x[i].end(); it++)
				S[i] = std::max(S[i], abs(*it));

			if (S[i] < eps)
			{
				std::cout << "Error: Matrix has a null column/row." << std::endl;
				return;
			}
			
			for (auto it = m_x[i].begin(); it != m_x[i].end(); it++)
				*it /= S[i];
		}
	}
	
	for (i = 0; i < ncols; i++)
		S[i]  = 1.0/S[i];
}

#endif