#ifndef _LILC_MATRIX_ILDL_H_
#define _LILC_MATRIX_ILDL_H_


using std::endl;
using std::cout;
using std::abs;

template <class el_type>
void lilc_matrix<el_type> :: ildl(lilc_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol)
{

	//----------------- initialize temporary variables --------------------//
	int lfil = 2*fill_factor*nnz()/n_cols(); //roughly a factor of 2 since only lower tri. of A is stored
	const el_type alpha = (1+sqrt(17))/8;  //for use in pivoting.
	el_type w1, wr, d1, dr(-1);
	el_type det_D, D_inv11, D_inv22, D_inv12;
	el_type l_11, l_12;

	const int ncols = n_cols(); //number of cols in A.

	vector<bool> in_set(ncols, false); //bitset used for unsorted merges
	swap_struct<el_type> s;	//struct containing temp vars used in pivoting.
	
	elt_vector_type work(ncols, 0), temp(ncols, 0); ////work vector for the current column
	idx_vector_type curr_nnzs, temp_nnzs;  //non-zeros on current col.
	curr_nnzs.reserve(ncols); //reserves space for worse case (entire col is non-zero)

	int count = 0; //the current non-zero in L.
	int i, j, k, r, offset, col_size, col_size2(-1);
	bool size_two_piv = false;	//boolean indicating if the pivot is 2x2 or 1x1

	//--------------- allocate memory for L and D ------------------//
	L.resize(ncols, ncols);
	L.list.resize(ncols ); //allocate a vector of size n for Llist.
	D.resize(ncols );
	

	int debug = -1;
	
	//------------------- main loop: factoring begins -------------------------//
	for (k = 0; k < ncols; k++) {

		//zero out work vector
		curr_nnzs.clear();

		if (k == debug) { 
			cout << "m_idx[k]: " << m_idx[k] << endl;
			cout << "m_x[k]: " << m_x[k] << endl;
		}
		
		//assign nonzeros indices of A(k:n, k) to curr_nnzs
		curr_nnzs.assign (m_idx[k].begin(), m_idx[k].end());

		//assign nonzero values of A(k:n, k) to work
		for (j = 0; j < (int) curr_nnzs.size(); j++) {
			work[curr_nnzs[j]] = m_x[k][j];
		}

		//--------------begin pivoting--------------//

		//do delayed updates on current column. work = Sum_{i=0}^{k-1} L(k,i) * L(k:n, i).
		update(k, work, curr_nnzs, L, D, in_set);
		
		//store diagonal element in d1. set diagonal element in work vector to 0
		//since we want to find the maximum off-diagonal element.
		d1 = work[k];
		work[k] = 0;

		//find maximum element in work and store its index in r.
		w1 = max(work, curr_nnzs, r);

		if (k == debug) { 
			cout << "curr_nnzs: " << curr_nnzs << endl;
			for (i = 0; i < (int) curr_nnzs.size(); i++) {
				cout << work[curr_nnzs[i]] << " ";
			}
			cout << endl;
			cout << "w1: " << w1 << endl;
		}
		
		//bunch-kaufman partial pivoting is used below. for a more detailed reference,
		//refer to "Accuracy and Stability of Numerical Algorithms." by Higham (2002).
		//------------------- begin bunch-kaufman pivoting ------------------//
		if (w1 < eps) {
			//case 0: do nothing. pivot is k.
			if (k == debug) {
				cout << "case 0" << endl;
			}
		} else if ( (alpha * w1 - abs(d1)) < eps  ) {
			//case 1: do nothing. pivot is k.
			if (k == debug) {
				cout << "case 1" << endl;
			}
		} else {
			//zero out temp vector.
			temp_nnzs.clear();

			//ensure invariant 1.
			ensure_invariant(r, k, list[r], first[r], true);
			offset = first[r];

			for (j = offset; j < (int) list[r].size(); j++) {
				temp_nnzs.push_back(list[r][j]);
				temp[list[r][j]] = coeff(r, list[r][j]);
			}

			unordered_inplace_union(temp_nnzs, m_idx[r].begin(), m_idx[r].end(), in_set);

			for (j = 0; j < (int) m_idx[r].size(); j++) {
				temp[m_idx[r][j]] = m_x[r][j];
			}

			update(r, temp, temp_nnzs, L, D, in_set);

			
			dr = temp[r];
			temp[r] = 0;

			wr = max(temp, temp_nnzs, j);

			if ((alpha*w1*w1 - abs(d1)*wr) < eps) {
				//case 2: do nothing. pivot is k.
				if (k == debug) {
					cout << "case 2" << endl;
					cout << "d1: " << d1 << " wr: " << wr << " w1: " << w1 << endl;
				}
				
				
			} else if ( (alpha * wr - abs(dr)) < eps) {
				//case 3: pivot is k with r: 1x1 pivot case.
				temp[r] = dr;
				work[k] = d1;

				if (k > r) {
					cout << "case 3 " << k << " " << r << endl;
					cout << "fault! " << k << " " << r << endl;
					return;
				}

				//--------pivot A and L ---------//
				pivot(s, in_set, L, k, r);

				//----------pivot rest ----------//
				std::swap(d1, dr);

				//permute perm
				std::swap(perm[k], perm[r]);

				work.swap(temp);	//swap work with temp.
				std::swap(work[k], work[r]);

				curr_nnzs.swap(temp_nnzs);	//swap curr_nnzs with temp_nnzs

				safe_swap(curr_nnzs, k, r);

				//--------end pivot rest---------//

			} else {
				//case 4: pivot is k+1 with r: 2x2 pivot case.
				if (k == debug) {
					cout << "case 4" << " " << k+1 << " " << r << endl;
					cout << "temp_nnzs: " << temp_nnzs << " " << r << endl;
					for (i = 0; i < (int) temp_nnzs.size(); i++) {
						cout << temp[temp_nnzs[i]] << " ";
					}
					cout << endl;
				}
				if (k >= r) {
					cout << "fault! " << k+1 << " " << r << endl;
					return;
				}

				advance_list(k);
				L.advance_first(k);

				temp[r] = dr;
				work[k] = d1;

				size_two_piv = true;

				if (k+1 != r) {
					pivot(s, in_set, L, k+1, r);


					//----------pivot rest ----------//

					//permute perm
					std::swap(perm[k+1], perm[r]);

					//swap two cols of L
					std::swap(work[k+1], work[r]);
					std::swap(temp[k+1], temp[r]);

					safe_swap(curr_nnzs, k+1, r);
					safe_swap(temp_nnzs, k+1, r);
				}
			}
		}
		//--------------end pivoting--------------//

		curr_nnzs.erase(std::remove(curr_nnzs.begin(), curr_nnzs.end(), k), curr_nnzs.end());

		//performs the dual dropping procedure.
		if (!size_two_piv) {
			drop_tol(work, curr_nnzs, lfil, tol);

		} else {
			temp_nnzs.erase(std::remove(temp_nnzs.begin(), temp_nnzs.end(), k), temp_nnzs.end());
			curr_nnzs.erase(std::remove(curr_nnzs.begin(), curr_nnzs.end(), k+1), curr_nnzs.end());
			temp_nnzs.erase(std::remove(temp_nnzs.begin(), temp_nnzs.end(), k+1), temp_nnzs.end());
			
			det_D = d1*dr - work[k+1]*work[k+1];
			if ( abs(det_D) < eps) det_D = 1e-3;  //statically pivot;
			D_inv11 = dr/det_D;
			D_inv22 = d1/det_D;
			D_inv12 = -work[k+1]/det_D;
			
			D.off_diagonal(k) = work[k+1];
			D[k+1] = dr;
			
			unordered_inplace_union(curr_nnzs, temp_nnzs.begin(), temp_nnzs.end(), in_set);

			
			for (auto it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) { //need -1 on col_size to remove offset from initializing col_size to 1
				l_11 = work[*it]*D_inv11 + temp[*it]*D_inv12;
				l_12 = work[*it]*D_inv12 + temp[*it]*D_inv22;
				
				work[*it] = l_11;
				temp[*it] = l_12;
			}
			
			temp_nnzs.assign(curr_nnzs.begin(), curr_nnzs.end());
			drop_tol(temp, temp_nnzs, lfil, tol);
			drop_tol(work, curr_nnzs, lfil, tol);
			

		}

		L.m_idx[k].resize(curr_nnzs.size()+1);
		L.m_x[k].resize(curr_nnzs.size()+1);
		
		D[k] = d1;
		
		L.m_x[k][0] = 1;
		L.m_idx[k][0] = k;
		count++;
		
		if (!size_two_piv) {
			if ( abs(D[k]) < eps) D[k] = 1e-3; //statically pivot
			i = 0;
			for (auto it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) { //need -1 on col_size to remove offset from initializing col_size to 1
				if ( abs(work[*it]) > eps) {
					L.m_idx[k][i+1] = *it; //row_idx of L is updated
					L.m_x[k][i+1] = work[*it]/D[k]; //work vector is scaled by D[k]

					L.list[*it].push_back(k); //update Llist
					count++;
					i++;
				}
			}
			
			col_size = 1 + i;
			//update lfirst
			L.advance_first(k);
			advance_list(k);
		} else {
			L.m_idx[k+1].resize(temp_nnzs.size()+1);
			L.m_x[k+1].resize(temp_nnzs.size()+1);
			
			D.off_diagonal(k) = work[k+1];
			D[k+1] = dr;

			L.m_x[k+1][0] = 1;
			L.m_idx[k+1][0] = k+1;
			count++;

			i = 0;
			for (auto it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) {
				if ( abs(work[*it]) > eps) {
					L.m_x[k][i+1] = work[*it]; //row_idx of L is updated
					L.m_idx[k][i+1] = *it; //work vector is scaled by D[k]
					L.list[*it].push_back(k); //update Llist
					count++;
					i++;
				}
				
			}
			
			j = 0;
			for (auto it = temp_nnzs.begin(); it != temp_nnzs.end(); it++) {
				if ( abs(temp[*it]) > eps) {
					L.m_x[k+1][j+1] = temp[*it]; //row_idx of L is updated
					L.m_idx[k+1][j+1] = *it; //work vector is scaled by D[k]
					L.list[*it].push_back(k+1); //update Llist
					count++;
					j++;
				}
				
			}

			col_size = 1 + i;
			col_size2 = 1 + j;

			//update lfirst
			L.advance_first(k+1);
			advance_list(k+1);
			
			
			// for (auto it = temp_nnzs.begin(); it != temp_nnzs.end(); it++) {
				// temp[*it] = 0;
			// }
			//}
		}
		
		//cleanup
		work[k] = 0;
		temp[k] = 0;
		
		if (k + 1 < ncols) {
			temp[k+1] = 0;
			work[k+1] = 0;
		}
		
		for (auto it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) {
			work[*it] = 0;
		}
		curr_nnzs.clear();
		
		for (auto it = temp_nnzs.begin(); it != temp_nnzs.end(); it++) {
			temp[*it] = 0;
		}
		temp_nnzs.clear();
		
		L.m_x[k].resize(col_size);
		L.m_idx[k].resize(col_size);

		if (k == debug) {
			cout << "part of two piv? " << size_two_piv << endl;
			cout << "lfil: " << lfil << endl;
			cout << "D[k]: " << D[k] << endl;
			cout << "D(k,k+1): " << D.off_diagonal(k) << endl;
			cout << "D[k+1]: " << D[k+1] << endl;
			
			cout << "curr_nnzs: " << curr_nnzs << " " << r << endl;
			for (i = 0; i < col_size; i++) {
				cout << work[curr_nnzs[i]]/D[k] << " ";
			}
			cout << endl;
		}
		
		if (size_two_piv) {
			L.m_x[k+1].resize(col_size2);
			L.m_idx[k+1].resize(col_size2);
			k++;
			
			size_two_piv = false;
		}
		
		// std::cout << k << std::endl;
		// std::string s;
		// std::cin >> s;
		// if (s == "l") {
		// std::cout << L << std::endl;
		// } else if (s == "a") {
		// std::cout << to_string() << std::endl;
		// } else if (s == "w") {
			// cout << work << endl;
		// } else if (s == "t") {
			// cout << temp << endl;
		// }
		
	}

	L.nnz_count = count;

}

#endif
