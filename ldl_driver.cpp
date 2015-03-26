#include "source/solver.h"

#include <iostream>
#include <cassert>
#include <cstring>
#include "include/gflags/gflags.h"

/*!	\mainpage \authors <a href="https://www.cs.ubc.ca/~greif/">Chen Greif</a>, Shiwen He, <a href="https://www.cs.ubc.ca/~inutard/">Paul Liu</a>
*
*	\tableofcontents
*	\section intro_sec Introduction
*
* 	
	\b sym-ildl is a C++ package for producing fast incomplete factorizations of symmetric indefinite matrices. Given an \f$n\times n\f$ symmetric indefinite matrix \f$\mathbf{A}\f$, this package produces an incomplete \f$\mathbf{LDL^{T}}\f$ factorization. Prior to factorization, this package first scales the matrix to be equilibriated in the max-norm [2], and then preorders the matrix using either the Reverse Cuthill-McKee (RCM) algorithm or the Approximate Minimum Degree algorithm (AMD) [1]. To maintain stability, the user can use Bunch-Kaufman or rook partial pivoting during the factorization process. The factorization produced is of the form 
	\f[
		\mathbf{P^{T}SASP=LDL^{T}}.
	\f]
	where \f$\mathbf{P}\f$ is a permutation matrix, \f$\mathbf{S}\f$ a scaling matrix, and \f$\mathbf{L}\f$ and \f$\mathbf{D}\f$ are the unit lower triangular and diagonal factors respectively. 
	
	This package is based on and extends an incomplete factorization approach proposed by Li and Saad [3] (which itself builds on Li, Saad, and Chow [4]).
	
	
*	\section quick_start Quick Start
*
	To begin using the package, first download the files hosted at <a href="https://github.com/inutard/matrix-factor">https://github.com/inutard/matrix-factor</a>. The package works under most Unix distributions as well as Cygwin under Windows. The default compiler used is \c g++, simply type \c make at the command line to compile the entire package. In addition to \subpage ldl_driver "usage as a standalone program", the package also has a \subpage matlab_mex "Matlab interface".

	\subsection ldl_driver Using the package as a standalone program
	The compiled program \c ldl_driver takes in (through the command line) three parameters as well as four optional ones.
	
	The format of execution is: 
	\code 
		./ldl_driver -filename=[matrix-name.mtx] -fill=[fill_factor] -tol=[drop_tol] -pp_tol=[pp_tol] -reordering=[amd/rcm/none] -save=[y/n] -display=[y/n]
	\endcode
	
	The parameters above can be given in any order, and will use a default value when not specified.
	
	A description of each of these parameters can be accessed by typing 
	\code 
		./ldl_driver --help 
	\endcode
	
	For convenience, the parameters are listed below:
	\param filename The filename of the matrix to be loaded. Several test matrices exist in the test_matrices folder. All matrices loaded are required to be in matrix market (.mtx) form.
	
	\param fill Controls memory usage. Each column is guaranteed to have fewer than \f$fill\cdot nnz(\mathbf{A})/n\f$ elements. When this argument is not given, the default value for \c fill is <c>1.0</c>.
	
	\param tol Controls agressiveness of dropping. In each column k, elements less than \f$tol \cdot \left|\left|\mathbf{L}_{k+1:n,k}\right|\right|_1\f$ are dropped. The default value for \c tol is <c>0.001</c>.
	
	\param pp_tol A parameter to aggressiveness of Bunch-Kaufman pivoting (BKP). When pp_tol >= 1, full BKP is used. When pp_tol is 0, there is no partial pivoting. Values between 0 and 1 varies the number of pivots of BKP makes.

	\param reordering Determines what sort of preordering will be used on the matrix. Choices are 'amd', 'rcm', and 'none'. The default is 'amd'.
	
	\param save Indicates whether the output matrices should be saved. \c y indicates yes, \c n indicates no. The default flag is \c y. All matrices are saved in matrix market (.mtx) form. The matrices are saved into an external folder named \c output_matrices. There are five saved files: <c>outA.mtx, outL.mtx, outD.mtx, outS.mtx</c>, and \c outP.mtx. \c outB.mtx is the matrix \f$\mathbf{B=P^{T}SASP}\f$. The rest of the outputs should be clear from the description above.
	
	\param display Indicates whether the output matrices should be displayed to the command line. \c y indicates yes, \c n indicates no. The default flag is \c y.	
	
	Typically, the \c pp_tol and \c reordering parameters are best left to the default options.
	
	\par Examples:
	Suppose we wish to factor the \c aug3dcqp matrix stored in <c>test_matrices/aug3dcqp.mtx</c>. Using the parameters described above, the execution of the program may go something like this:	
	\code
		./ldl_driver -filename=test_matrices/aug3dcqp.mtx -fill=1.0 tol=0.001 -save=y -display=n
		
		Load succeeded. File test_matrices/aug3dcqp.mtx was loaded.
		A is 35543 by 35543 with 128115 non-zeros.
		Equilibration:  0.047 seconds.
		AMD:            0.047 seconds.
		Permutation:    0.047 seconds.
		Factorization:  0.109 seconds.
		Total time:     0.250 seconds.
		L is 35543 by 35543 with 108794 non-zeros.
		Saving matrices...
		Save complete.
		Factorization Complete. All output written to /output_matrices directory.
	\endcode	
	The code above factors the \c aug3dcqp.mtx matrix (<c>lfil=1.0, tol=0.001</c>) from the \c test_matrices folder and saves the outputs. The time it took to pre-order and equilibriate the matrix (0.047s) as well as the actual factorization (0.109s) are also given.
	
	\par
	The program may also run without the last 4 arguments:	
	\code
		./ldl_driver -filename=test_matrices/aug3dcqp.mtx -fill=1.0 -tol=0.001
		
		Load succeeded. File test_matrices/aug3dcqp.mtx was loaded.
		A is 35543 by ...
	\endcode
	This code does the exact same thing as the code in the previous example, except this time we take advantage of the fact that \c save defaults to \c y and \c display to \c n.
	
	\par
	Finally, we may use all optional arguments:
	\code
		./ldl_driver -filename=test_matrices/aug3dcqp.mtx
		
		Load succeeded. File test_matrices/aug3dcqp.mtx was loaded.
		A is 35543 by ...
	\endcode
	The code above would use the default arguments <c>-fill=1.0 -tol=0.001 -pp_tol=1.0 -reordering=amd -save=y -display=n</c>.
	
	\subsection matlab_mex Using sym-ildl within Matlab
	If everything is compiled correctly, simply open Matlab in the package directory. The \c startup.m script adds all necessary paths to Matlab upon initiation. The program can now be called by its function handle, \c ildl.
	
	\c ildl takes in five arguments, four of them being optional. A full description of the parameters can be displayed by typing 
	\code
		help ildl
	\endcode
	
	For convenience, the parameters are listed below:
	\param A The matrix to be factored.
	
	\param fill Controls memory usage. Each column is guaranteed to have fewer than \f$fill\cdot nnz(\mathbf{A})/n\f$ elements. When this argument is not given, the default value for \c fill is <c>1.0</c>.
	
	\param tol Controls agressiveness of dropping. In each column k, elements less than \f$tol \cdot \left|\left|\mathbf{L}_{k+1:n,k}\right|\right|_1\f$ are dropped. The default value for \c tol is <c>0.001</c>.
	
	\param pp_tol A parameter to aggressiveness of Bunch-Kaufman pivoting (BKP). When pp_tol >= 1, full BKP is used. When pp_tol is 0, there is no partial pivoting. Values between 0 and 1 varies the number of pivots of BKP makes.

	\param reordering Determines what sort of preordering will be used on the matrix. Choices are 'amd', 'rcm', and 'none'. The default is 'amd'.
	
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
		>> B = sparse(gallery('uniformdata',100,0));
		>> A = [speye(100) B; B' sparse(100, 100)];
	\endcode
	The \c A generated is a special type of matrix called a KKT matrix. These matrices are indefinite and arise often in optimzation problems. Note that A must be a Matlab \b sparse matrix.
	
	\par
	To factor the matrix, we supply \c ildl with the parameters described above:
	\code
		>> [L, D, p, S, B] = ildl(A, 1.0, 0.001);
		Equilibration:	0.001 seconds.
		AMD:		0.001 seconds.
		Permutation:	0.000 seconds.
		Factorization:	0.022 seconds.
		Total time:	0.024 seconds.
		L is 200 by 200 with 14388 non-zeros.
	\endcode
	As we can see above, \c ildl will supply some timing information to the console when used. The reordering time is the time taken to equilibriate and preorder the matrix. The factorization time is the time it took to factor and pivot the matrix with partial pivoting.
	
	\par
	We may also take advantage of the optional parameters and simply feed \c ildl only one parameter:
	\code
		>> [L, D, p, S, B] = ildl(A);
		Equilibration:	0.001 seconds.
		AMD:		0.001 seconds.
		...
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

DEFINE_string(filename, "", "The filename of the matrix to be factored"
							"(in matrix-market format).");
							
DEFINE_double(fill, 1.0, "A parameter to control memory usage. Each column is guaranteed"
						 "to have fewer than fill*nnz(A) elements.");
						 
DEFINE_double(tol, 0.001, "A parameter to control agressiveness of dropping. In each column k,"
						  "elements less than tol*||L(k+1:n,k)|| (1-norm) are dropped.");
						  
DEFINE_double(pp_tol, 1.0, "A parameter to aggressiveness of Bunch-Kaufman pivoting (BKP). "
						   "When pp_tol >= 1, full BKP is used. When pp_tol is 0, no BKP"
						   "is used. Values between 0 and 1 varies the aggressivness of"
						   "BKP in a continuous manner.");

DEFINE_string(pivot, "rook", "Determines what kind of pivoting algorithm will be used"
								 " during the factorization. Choices are 'rook' and 'bunch'.");
                                 
DEFINE_string(reordering, "amd", "Determines what sort of preordering will be used"
								 " on the matrix. Choices are 'amd', 'rcm', and 'none'.");
								 
DEFINE_bool(equil, true, "Decides if the matrix should be equilibriated before factoring is done. "
						 "If yes, matrix is equilibrated with Bunch's algorithm in the max norm.");

DEFINE_bool(inplace, false, "Decides if the matrix should be factored in place (faster and saves memory, "
						    "at the cost of not being able to use the built-in solver).");
						 
DEFINE_bool(save, true, "If yes, saves the factors (in matrix-market format) into a folder "
						 "called output_matrices/ in the same directory as ldl_driver.");

DEFINE_bool(display, false, "If yes, outputs a human readable version of the factors onto"
							" standard out. Generates a large amount of output if the "
							"matrix is big.");
							
DEFINE_int32(minres_iters, -1, "If >= 0 and supplied with a right hand side, SYM-ILDL will "
							   "attempt to use the preconditioner generated with MINRES to "
							   "solve the system.");
							   
DEFINE_string(rhs_file, "", "The filename of the right hand side (in matrix-market format).");

int main(int argc, char* argv[])
{
	std::string usage("Performs an incomplete LDL factorization of a given matrix.\n"
				      "Sample usage:\n"
					  "\t./ldl_driver -filename=test_matrices/testmat1.mtx "
									 "-fill=2.0 -display=true -save=false\n"
					  "Additionally, these flags can be loaded from a single file "
					  "with the option -flagfile=[filename].");
				 
	google::SetUsageMessage(usage);
	google::ParseCommandLineFlags(&argc, &argv, true);
	
	if (FLAGS_filename.empty()) {
		std::cerr << "No file specified! Type ./ldl_driver --help for a description of the program parameters." << std::endl;
		return 0;
	}
	
	solver<double> solv;
	solv.load(FLAGS_filename);
	
	//default reordering scheme is AMD
	solv.set_reorder_scheme(FLAGS_reordering.c_str());
	
	//default is equil on
	solv.set_equil(FLAGS_equil); 
	
	if (FLAGS_minres_iters > 0) {
        vector<double> rhs;
        if (!FLAGS_rhs_file.empty()) {
            read_vector(rhs, FLAGS_rhs_file);
        } else {
            // for testing purposes only
            rhs.resize(solv.A.n_cols(), 1);
        }
		
		if (rhs.size() != solv.A.n_cols()) {
			std::cout << "The right hand side dimensions do not match the dimensions of A." << std::endl;
			return 1;
		}
		solv.set_rhs(rhs);
	}
	
    solv.set_pivot(FLAGS_pivot.c_str());
    solv.set_inplace(FLAGS_inplace);
	solv.solve(FLAGS_fill, FLAGS_tol, FLAGS_pp_tol, FLAGS_minres_iters);
	
	if (FLAGS_save) {
		solv.save();
	}
	
	if (FLAGS_display) {
		solv.display();
		std::cout << endl;
	}
    
	std::cout << "Factorization Complete. All output written to /output_matrices directory." << std::endl;

	return 0;
}
