//-*- mode: c++ -*-
#ifndef _SOLVER_MINRES_H_
#define _SOLVER_MINRES_H_

#include <string>
#include <sstream>

template <class el_type>
void solver<el_type> :: minres(int max_iter, double shift) {
	int n = A.n_rows();
	sol_vec.resize(n, 0);
	
	// ---------- set initial values and allocate memory for variables ---------//
	el_type alpha[2], beta[2]; // the last two entries of the T matrix
	el_type res[2]; // the last two residuals
	
	// the last two vectors of the lanczos iteration
	vector< vector<el_type> > v(vector<el_type>(n), 2);
	
	double norm_A = 0; // matrix norm estimate
	double cond_A = 1;	// condition number estimate
	double c = 1, s = 0; // givens rotation elements
	
	// temporary variables to store the corner of the matrix we're factoring
	double gamma_min = 1;
	el_type delta1[2], delta2[2], ep[2], gamma1[2];
	
	// temporary vectors for lanczos calcluations
	vector<el_type> pk(n);
	
	// step size in the current search direction (xk = x_{k-1} + tau*dk)
	double tau = 0;
	
	// the last 3 search directions
	vector< vector<el_type> > d(vector<el_type>(n), 3);
	
	// set up initial values for variables above
	double eps = A::eps;
	beta[0] = norm(rhs, 2);
	
	// v[0] = rhs/beta[0]
	for (int i = 0; i < n; i++) {
		v[0][i] = rhs[i]/beta[0];
	}
	
	res[0] = beta[0];
	tau = beta[0];
	
	auto sign = [&] { return (abs(x) < eps ? 0 : x/abs(x)); };
	
	// -------------- begin minres iterations --------------//
	int k = 0; // iteration number
	while (res[cur] < eps && k < max_iter) {
		int cur = k%2, nxt = (k+1)%2;
		// ---------- begin lanczos step ----------//
		//pk = (A - shift*I) * v[cur]
		A.multiply(v[cur], pk);
		for (int i = 0; i < n; i++) {
			pk[i] -= shift;
		}
		
		// alpha = v[cur]' * pk
		alpha[cur] = dot_product(v[cur], pk);
		
		// pk = pk - alpha*v[cur]
		vector_sum(1, pk, -alpha[cur], v[cur], pk);
		
		// v[nxt] = pk - beta[cur]*v[nxt]
		vector_sum(1, pk, -beta[cur], v[nxt], v[nxt]);
		beta[nxt] = norm(v[nxt], 2);
		
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
		ep[nxt] = s*beta[next];
		delta1[nxt] = -c*beta[nxt];
		
		// ---------- begin givens rotation ----------//
		double a = gamma1[cur], b = beta[nxt];
		if (abs(b) < eps) {
			s = 0; 
			gamma2[nxt] = abs(a);
			if (abs(a) < eps) {
				c = 1;
			} else {
				c = sign(a);
			}
		} else if (abs(a) < eps) {
			c = 0;
			s = sign(b);
			gamma2[nxt] = abs(b);
		} else if (abs(b) > abs(a)) {
			double t = a/b;
			s = sign(b)/sqrt(1+t*t);
			c = s*t;
			gamma2[nxt] = b/s;
		} else { //abs(a) >= abs(b)
			double t = b/a;
			c = sign(a)/sqrt(1+t*t);
			s = c*t;
			gamma2[nxt] = a/c;
		}
		// ---------- end givens rotation ----------//
		
		// update residual norms and estimate for matrix norm
		tau = c*res[cur];
		res[nxt] = s*res[cur];
		
		double tnorm = sqrt(alpha[cur]*alpha[cur] + beta[nxt]*beta[nxt]);
		norm_A = max(norm_A, tnorm);
		
		// ------ update solution and matrix condition number ------ //
		if (abs(gamma2[cur]) > eps) {
			// compute new search direction
			// d[(k+2)%3] = (v[cur] - delta2[cur]*d[(k+1)%3] - ep[cur]*d[k%3])/gamma2[cur];
			vector_sum(1, v[cur], -delta2[cur], d[(k+1)%3], d[(k+2)%3]);
			vector_sum(1, d[(k+2)%3], -ep[cur], d[k%3], d[(k+2)%3]);
			
			for (int i = 0; i < n; i++) {
				d[(k+2)%3][i] /= gamma2[cur];
			}
			
			//sol = sol + tau*d[(k+2)%3]
			vector_sum(1, sol_vec, tau, d[(k+2)%3]);
			gamma_min = min(gamma_min, gamma2[cur]);
			cond_A = norm_A/gamma_min;
		}
		
		k++;
		
		// ------------- end update ------------- //
	}
	
	printf("The estimated condition number of the matrix is %e.\n", cond_A);
	printf("Minres took %i iterations and got down to %e residual.\n", k, res[k%2]);
	
}

#endif // _SOLVER_MINRES_H_
