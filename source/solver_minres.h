//-*- mode: c++ -*-
#ifndef _SOLVER_MINRES_H_
#define _SOLVER_MINRES_H_

#include <algorithm>

template<class el_type, class mat_type >
void solver<el_type, mat_type> :: minres(int max_iter, double stop_tol, double shift) {
	int n = A.n_rows();
	sol_vec.resize(n, 0);
	
	// ---------- set initial values and allocate memory for variables ---------//
	el_type alpha[2], beta[2]; // the last two entries of the T matrix
	el_type res[2]; // the last two residuals
	
	// the last two vectors of the lanczos iteration
	vector< vector<el_type> > v(2, vector<el_type>(n));
	
	double norm_A = 0; // matrix norm estimate
	double cond_A = 1;	// condition number estimate
	double c = -1, s = 0; // givens rotation elements
	
	// temporary variables to store the corner of the matrix we're factoring
	double gamma_min = 1;
	el_type delta1[2], delta2[2], ep[2], gamma1[2], gamma2[2];

	// temporary vectors for lanczos calcluations
	vector<el_type> pk(n), tk(n);
	
	// step size in the current search direction (xk = x_{k-1} + tau*dk)
	double tau = 0;
	
	// the last 3 search directions
	vector< vector<el_type> > d(2, vector<el_type>(n));
	
	// set up initial values for variables above
	double eps = A.eps;
	beta[0] = norm(rhs, 2.0);
	
	double norm_rhs = beta[0];
	
	// v[0] = rhs/beta[0]
	for (int i = 0; i < n; i++) {
		v[0][i] = rhs[i]/beta[0];
	}
	
	res[0] = beta[0];
	tau = beta[0];
	
	auto sign = [&](double x) { return (abs(x) < eps ? 0 : x/abs(x)); };

	// -------------- begin minres iterations --------------//
	int k = 0; // iteration number
	while (res[k%2]/norm_rhs > stop_tol && k < max_iter) {
		int cur = k%2, nxt = (k+1)%2;
		// ---------- begin lanczos step ----------//
		//pk = (M^(-1) A M^(-t) - shift*I) * v[cur], where M = L|D|^(1/2) where |D|^(1/2) = Q|V|^(1/2)
		//we do this in steps. first, tk = L^(-t) * |D|^(-t/2) 
		D.sqrt_solve(v[cur], pk, true);
		L.forwardsolve(pk, tk);
		
		//pk = A*tk
		A.multiply(tk, pk);

		//pk = |D|^(-1/2) L^(-1) pk. after this step, pk = M^(-1) A M^(-t) v[cur]
		L.backsolve(pk, tk);
		D.sqrt_solve(tk, pk, false);
		
		//finally, pk = pk - shift*I * v[cur];
		for (int i = 0; i < n; i++) {
			pk[i] -= shift * v[cur][i];
		}
		
		// alpha = v[cur]' * pk
		alpha[cur] = dot_product(v[cur], pk);
		
		// pk = pk - alpha*v[cur]
		vector_sum(1, pk, -alpha[cur], v[cur], pk);
		
		// v[nxt] =  pk - beta[cur]*v[nxt]
		vector_sum(1, pk, -beta[cur], v[nxt], v[nxt]);
		beta[nxt] = norm(v[nxt], 2.0);
		
		// scale v[nxt] if beta[nxt] is not zero
		if (abs(beta[nxt]) > eps) {
			for (int i = 0; i < n; i++) {
				v[nxt][i] /= beta[nxt];
			}
		}
		// ---------- end lanczos step ----------//
		
		// left orthogonlization on the middle two entries in the last column of Tk
		delta2[cur] = c*delta1[cur] + s*alpha[cur];
		gamma1[cur] = s*delta1[cur] - c*alpha[cur];
		
		// left orthogonalization to product first two entries of T_{k+1} and ep_{k+1}
		ep[nxt] = s*beta[nxt];
		delta1[nxt] = -c*beta[nxt];
		
		// ---------- begin givens rotation ----------//
		double a = gamma1[cur], b = beta[nxt];
		if (abs(b) < eps) {
			s = 0; 
			gamma2[cur] = abs(a);
			if (abs(a) < eps) {
				c = 1;
			} else {
				c = sign(a);
			}
		} else if (abs(a) < eps) {
			c = 0;
			s = sign(b);
			gamma2[cur] = abs(b);
		} else if (abs(b) > abs(a)) {
			double t = a/b;
			s = sign(b)/sqrt(1+t*t);
			c = s*t;
			gamma2[cur] = b/s;
		} else { //abs(a) >= abs(b)
			double t = b/a;
			c = sign(a)/sqrt(1+t*t);
			s = c*t;
			gamma2[cur] = a/c;
		}
		// ---------- end givens rotation ----------//
		
		// update residual norms and estimate for matrix norm
		tau = c*res[cur];
		res[nxt] = s*res[cur];
		
		if (k == 0) norm_A = sqrt(alpha[cur]*alpha[cur] + beta[nxt]*beta[nxt]);
		else {
			double tnorm = sqrt(alpha[cur]*alpha[cur] + beta[nxt]*beta[nxt] + beta[cur]*beta[cur]);
			norm_A = std::max(norm_A, tnorm);
		}
		
		// ------ update solution and matrix condition number ------ //
		if (abs(gamma2[cur]) > eps) {
			// compute new search direction
			// d[nxt] = (v[cur] - delta2[cur]*d[cur] - ep[cur]*d[nxt])/gamma2[cur];
			vector_sum(1, v[cur], -ep[cur], d[nxt], d[nxt]);
			vector_sum(1, d[nxt], -delta2[cur], d[cur], d[nxt]);
			
			for (int i = 0; i < n; i++) {
				d[nxt][i] /= gamma2[cur];
			}
			
			//sol = sol + tau*d[nxt]
			vector_sum(1, sol_vec, tau, d[nxt], sol_vec);
			gamma_min = std::min(gamma_min, gamma2[cur]);
			cond_A = norm_A/gamma_min;
		}
		
		k++;
		
		// ------------- end update ------------- //
		
		//cout << "current residual " << res[cur]/norm_rhs << endl;
	}
	
	printf("The estimated condition number of the matrix is %e.\n", cond_A);
	printf("MINRES took %i iterations and got down to relative residual %e.\n", k, res[k%2]/norm_rhs);
	return;
}

#endif // _SOLVER_MINRES_H_
