#include <iostream>
#include <cassert>
#include <string.h>
#include "source/solver.h"

using namespace std;

/*! \mainpage Main Page
*
* \section intro_sec Introduction
*
* 	\b matrix-factor is a C++ package for producing fast incomplete factorizations of symmetric indefinite matrices. Given an \f$n\times n\f$ symmetric indefinite matrix \f$\mathbf{A}\f$, this package produces an incomplete \f$\mathbf{LDL^{T}}\f$ factorization. Prior to factorization, this package first scales the matrix to be equilibriated in the max-norm, and then preorders the matrix using the Reverse Cuthill-McKee algorithm. To maintain stability, we use Bunch-Kaufman partial pivoting during the factorization process. The factorization produced is of the form 
	\f[
		\mathbf{P^{T}SASP=LDL^{T}}.
	\f]
	where \f$\mathbf{P}\f$ is a permutation matrix, \f$\mathbf{S}\f$ a scaling matrix, and \f$\mathbf{L}\f$ and \f$\mathbf{D}\f$ are the unit lower triangular and diagonal factors respectively. 
	
*	\section quick_start Quick Start
*
*	To begin using the package, first download the files hosted at <a href="https://github.com/inutard/matrix-factor">https://github.com/inutard/matrix-factor</a>. The package works under most Unix distributions as well as Cygwin under Windows. The default compiler used is \c gcc, simply type \c make at the command line to compile.

	The compiled program \c ldl_driver takes in (through the command line) three parameters as well as two optional ones.
	
	The format of execution is: 
	\code 
	./ldl_driver [fill_factor] [tol] [in.mtx] [save] [display]
	\endcode
	
	The parameters are listed below:
	\param fill_factor A parameter to control memory usage. Each column is guaranteed to have fewer than \f$fill\_factor\cdot nnz(\mathbf{A})/n\f$ elements.
	
	\param tol A parameter to control agressiveness of dropping. In each column k, elements less than \f$tol \cdot \left|\left|\mathbf{L}_{k+1:n,k}\right|\right|_1\f$ are dropped.
	
	\param in.mtx The filename of the matrix to be loaded. Several test matrices exist in the test_matrices folder. All matrices loaded are required to be in matrix market (.mtx) form.
	
	\param save A flag indicating whether the output matrices should be saved. \c -y indicates yes, \c -n indicates no. When this argument is not used, the default flag is \c -y. All matrices are saved in matrix market (.mtx) form. The matrices are saved into an external folder named \c output_matrices. There are five saved files: <c>outA.mtx, outL.mtx, outD.mtx, outS.mtx</c>, and \c outP.mtx. \c outB.mtx is the matrix \f$\mathbf{B=P^{T}SASP}\f$. The rest of the outputs should be clear from the description above.
	
	\param display A flag indicating whether the output matrices should be displayed to the command line. \c -y indicates yes, \c -n indicates no. When this argument is not used, the default flag is \c -n.
	
	\par Examples:
	Using the parameters described above, the execution of the program may go something like this:	
	\code
	./ldl_driver 1.0 0.001 test_matrices/testmat1.mtx -y -y
	\endcode	
	The code above factors the \c testmat1.mtx matrix (<c>lfil=1.0, tol=0.001</c>) from the \c test_matrices folder, displays the full factorization (L and D) to terminal, and saves the outputs. 
	The program may also run without the last two arguments:	
	\code
	./ldl_driver 1.0 0.001 test_matrices/testmat1.mtx
	\endcode
	This code uses the default flags <c>-y -n</c> for the last two arguments, resulting in the outputs being saved, but not displayed to the terminal.

*	\section refs References
	-#	J. A. George and J. W-H. Liu, "Computer Solution of Large Sparse Positive Definite Systems", 1981.
	-#	J. R. Bunch, "Equilibration of Symmetric Matrices in the Max-Norm", 1971.
	-#	N. Li and Y. Saad, "Crout versions of the ILU factorization with pivoting for sparse symmetric matrices", 2005.
	-#	N. Li, Y. Saad, and E. Chow, "Crout versions of ILU for general sparse matrices", 2003.

*	
*/

int main(int argc, char* argv[]) {

	if (argc < 4 || argc > 6) {
		std::cout << "Too many or too few arguments." << std::endl
		<< "Program usage: ./ldl_driver [lfil] [tol] [in.mtx] [-save] [-display]" << std::endl;
		return 0;
	}

	solver<double> solv;
	solv.load(argv[3]);
	solv.solve(atof(argv[1]), atof(argv[2]));
	
	if (argc > 4) {
		if (strcmp(argv[4], "-y") == 0) {
			solv.save();
		}
		
		if (argc > 5 && strcmp(argv[5], "-y") == 0) {
			solv.display();
		}
		cout << "Done." << endl;
	} else {
		solv.save();
	}

	return 0;
}
