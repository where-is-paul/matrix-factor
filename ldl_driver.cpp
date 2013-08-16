#include <iostream>
#include <cassert>
#include <string.h>
#include "source/solver.h"

using namespace std;

/*!	\mainpage Main Page
*
*	
*	\section intro_sec Introduction
*
* 	
	\b sym-ildl is a C++ package for producing fast incomplete factorizations of symmetric indefinite matrices. Given an \f$n\times n\f$ symmetric indefinite matrix \f$\mathbf{A}\f$, this package produces an incomplete \f$\mathbf{LDL^{T}}\f$ factorization. Prior to factorization, this package first scales the matrix to be equilibriated in the max-norm, and then preorders the matrix using the Reverse Cuthill-McKee (RCM) algorithm. To maintain stability, we use Bunch-Kaufman partial pivoting during the factorization process. The factorization produced is of the form 
	\f[
		\mathbf{P^{T}SASP=LDL^{T}}.
	\f]
	where \f$\mathbf{P}\f$ is a permutation matrix, \f$\mathbf{S}\f$ a scaling matrix, and \f$\mathbf{L}\f$ and \f$\mathbf{D}\f$ are the unit lower triangular and diagonal factors respectively. 
	
*	\section quick_start Quick Start
*
	To begin using the package, first download the files hosted at <a href="https://github.com/inutard/matrix-factor">https://github.com/inutard/matrix-factor</a>. The package works under most Unix distributions as well as Cygwin under Windows. The default compiler used is \c gcc, simply type \c make at the command line to compile the entire package. Currently this package has only been tested for \c gcc-4.6.1 and above. In addition to \subpage ldl_driver "usage as a standalone program", the package also has a \subpage matlab_mex "Matlab interface".

	\subsection ldl_driver Using the package as a standalone program
	The compiled program \c ldl_driver takes in (through the command line) three parameters as well as two optional ones.
	
	The format of execution is: 
	\code 
		./ldl_driver [in.mtx] [fill_factor] [tol] [save] [display]
	\endcode
	
	The parameters are listed below:
	\param in.mtx The filename of the matrix to be loaded. Several test matrices exist in the test_matrices folder. All matrices loaded are required to be in matrix market (.mtx) form.
	
	\param fill_factor A parameter to control memory usage. Each column is guaranteed to have fewer than \f$fill\_factor\cdot nnz(\mathbf{A})/n\f$ elements. When this argument is not given, the default value for \c fill_factor is <c>1.0</c>.
	
	\param tol A parameter to control agressiveness of dropping. In each column k, elements less than \f$tol \cdot \left|\left|\mathbf{L}_{k+1:n,k}\right|\right|_1\f$ are dropped. The default value for \c tol is <c>0.001</c>.
	
	\param save A flag indicating whether the output matrices should be saved. \c -y indicates yes, \c -n indicates no. The default flag is \c -y. All matrices are saved in matrix market (.mtx) form. The matrices are saved into an external folder named \c output_matrices. There are five saved files: <c>outA.mtx, outL.mtx, outD.mtx, outS.mtx</c>, and \c outP.mtx. \c outB.mtx is the matrix \f$\mathbf{B=P^{T}SASP}\f$. The rest of the outputs should be clear from the description above.
	
	\param display A flag indicating whether the output matrices should be displayed to the command line. \c -y indicates yes, \c -n indicates no. The default flag is \c -n.
	
	\par Examples:
	Suppose we wish to factor the \c aug3dcqp matrix stored in <c>test_matrices/aug3dcqp.mtx</c>. Using the parameters described above, the execution of the program may go something like this:	
	\code
		./ldl_driver test_matrices/aug3dcqp.mtx 1.0 0.001 -y -y
		
		Load succeeded. File test_matrices/aug3dcqp.mtx was loaded.
		A is 35543 by 35543 with 77829 non-zeros.
		The reordering took 0.108006 seconds.
		The factorization took 0.101006 seconds.
		L is 35543 by 35543 with 143512 non-zeros.
		Saving matrices...
		Save complete.
	\endcode	
	The code above factors the \c aug3dcqp.mtx matrix (<c>lfil=1.0, tol=0.001</c>) from the \c test_matrices folder, displays the full factorization (L and D) to terminal, and saves the outputs. The time it took to pre-order and equilibriate the matrix (0.108s) as well as the actual factorization (0.101s) are also given.
	
	\par
	The program may also run without the last two arguments:	
	\code
		./ldl_driver test_matrices/aug3dcqp.mtx 1.0 0.001
		
		Load succeeded. File test_matrices/aug3dcqp.mtx was loaded.
		A is 35543 by ...
	\endcode
	This code uses the default flags <c>-y -n</c> for the last two arguments, resulting in the outputs being saved, but not displayed to the terminal.
	
	\par
	Finally, we may use all optional arguments:
	\code
		./ldl_driver test_matrices/aug3dcqp.mtx
		
		Load succeeded. File test_matrices/aug3dcqp.mtx was loaded.
		A is 35543 by ...
	\endcode
	The code above would use the default arguments <c>1.0 0.001 -y -n</c>.
	
	\subsection matlab_mex Using the package within Matlab
	If everything is compiled correctly, simply open Matlab in the package directory. The \c startup.m script adds all necessary paths to Matlab upon initiation. The program can now be called by its function handle, \c ildl.
	
	\c ildl takes in three arguments, two of them being optional. The parameters are listed below:
	\param A The matrix to be factored.
	
	\param fill_factor A parameter to control memory usage. Each column is guaranteed to have fewer than \f$fill\_factor\cdot nnz(\mathbf{A})/n\f$ elements. When this argument is not given, the default value for \c fill_factor is <c>1.0</c>.
	
	\param tol A parameter to control agressiveness of dropping. In each column k, elements less than \f$tol \cdot \left|\left|\mathbf{L}_{k+1:n,k}\right|\right|_1\f$ are dropped. The default value for \c tol is <c>0.001</c>.
	
	As with the standalone executable, the function has five outputs: <c>L, D, p, S,</c> and \c B:
	\return \b L Unit lower triangular factor of \f$\mathbf{P^{T}SASP}\f$.
	\return \b D Block diagonal factor (consisting of 1x1 and 2x2 blocks) of \f$\mathbf{P^{T}SASP}\f$.
	\return \b p Permutation vector containing permutations done to \f$\mathbf{A}\f$.
	\return \b S Diagonal scaling matrix that equilibrations \f$\mathbf{A}\f$ in the max-norm.
	\return \b B Permuted and scaled matrix \f$\mathbf{B=P^{T}SASP}\f$ after factorization.
*
*	\par Examples:
	Before we begin, let's first generate some symmetric indefinite matrices:
	\code
		>> B = sparse(gallery('uniformdata',10,0));
		>> A = [speye(10) B; B' sparse(10, 10)];
	\endcode
	The \c A generated is a special type of matrix called a KKT matrix. These matrices are indefinite and arise often in optimzation problems. Note that A must be a Matlab \b sparse matrix.
	
	\par
	To factor the matrix, we supply \c ildl with the parameters described above:
	\code
		>> [L, D, p, S, B] = ildl(A, 1.0, 0.001);
		The reordering took 0.000109 seconds.
		The factorization took 0.000456 seconds.
		L is 20 by 20 with 165 non-zeros.
	\endcode
	As we can see above, \c ildl will supply some timing information to the console when used. The reordering time is the time taken to equilibriate and preorder the matrix. The factorization time is the time it took to factor and pivot the matrix with partial pivoting.
	
	\par
	We may also take advantage of the optional parameters and simply feed \c ildl only one parameter:
	\code
		>> [L, D, p, S, B] = ildl(A);
		The reordering took 0.000109 seconds.
		The factorization took 0.000456 seconds.
		L is 20 by 20 with 165 non-zeros.
	\endcode

*
*
*	\section refs References
	-#	J. A. George and J. W-H. Liu, <em>Computer Solution of Large Sparse Positive Definite Systems</em>, Prentice-Hall, 1981.
	-#	J. R. Bunch, <em>Equilibration of Symmetric Matrices in the Max-Norm</em>, JACM, 18 (1971), pp. 566-572.
	-#	N. Li and Y. Saad, <em>Crout versions of the ILU factorization with pivoting for sparse symmetric matrices</em>, ETNA, 20 (2006), pp. 75-85.
	-#	N. Li, Y. Saad, and E. Chow, <em>Crout versions of ILU for general sparse matrices</em>, SISC, 25 (2003), pp. 716-728.

*	
*/

int main(int argc, char* argv[]) {

	if (argc < 2 || argc > 6) {
		std::cout << "Too many or too few arguments." << std::endl
		<< "Program usage: ./ldl_driver [in.mtx] [fill_factor] [tol] [-save] [-display]" << std::endl;
		std::cout << "Sample usage: ./ldl_driver test_matrices/testmat1.mtx 1.0 0.001 -y -y" << endl;
		return 0;
	}

	solver<double> solv;
	solv.load(argv[1]);
	
	double fill = 1.0, tol = 0.001;
	if (argc > 2) fill = atof(argv[2]);
	if (argc > 3) tol = atof(argv[3]);
	
	solv.solve(fill, tol);
	
	if (argc > 4) {
		if (strcmp(argv[4], "-y") == 0) {
			solv.save();
		}
		
		if (argc > 5 && strcmp(argv[5], "-y") == 0) {
			solv.display();
			cout << endl;
		}
		cout << "Factorization Complete. All output written to /output_matrices directory." << endl;
	} else {
		solv.save();
	}

	return 0;
}
