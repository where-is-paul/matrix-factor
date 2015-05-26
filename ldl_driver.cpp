#include "source/solver.h"

#include <iostream>
#include <cassert>
#include <cstring>
#include "include/gflags/gflags.h"

DEFINE_string(filename, "", "The filename of the matrix to be factored"
							"(in matrix-market format).");
							
DEFINE_double(fill, 3.0, "A parameter to control memory usage. Each column is guaranteed"
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

DEFINE_bool(full, false, "Uses SYM-ILDL to do the full factorization of the matrix, hence "
						 "using it as a direct solver.");
                            
DEFINE_bool(save, true, "If yes, saves the factors (in matrix-market format) into a folder "
						 "called output_matrices/ in the same directory as ldl_driver.");

DEFINE_bool(display, false, "If yes, outputs a human readable version of the factors onto"
							" standard out. Generates a large amount of output if the "
							"matrix is big.");
							
DEFINE_int32(minres_iters, -1, "If >= 0 and supplied with a right hand side, SYM-ILDL will "
							   "attempt to use the preconditioner generated with MINRES to "
							   "solve the system. The solution will be written to "
                               "output_matrices/ in matrix-market format.");

DEFINE_double(minres_tol, 1e-6, "A tolerance for MINRES. When the iterate x satisfies ||Ax-b||/||b|| < minres_tol, MINRES is terminated.");
                           
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
	
	if (FLAGS_minres_iters > 0 || FLAGS_full || !FLAGS_rhs_file.empty()) {
        if (!FLAGS_full && FLAGS_minres_iters <= 0) {
            printf("Using MINRES (200 max iterations) as default solver since RHS was loaded.\n");
            FLAGS_minres_iters = 200;
        }
        
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
	solv.set_full_solve(FLAGS_full);
    solv.solve(FLAGS_fill, FLAGS_tol, FLAGS_pp_tol, FLAGS_minres_iters, FLAGS_minres_tol);
	
	if (FLAGS_save) {
		solv.save();
	}
	
	if (FLAGS_display) {
		solv.display();
		std::cout << endl;
	}
    
	std::cout << "Factorization Complete. ";
    if (FLAGS_save) {
        std::cout << "All output written to /output_matrices directory.";
    }
    std::cout << std::endl;

	return 0;
}
