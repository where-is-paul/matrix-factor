#include "skew_symmetry_solver.h"

#include <iostream>
#include <cassert>
#include <cstring>
#include "gflags.h"

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
						   
DEFINE_string(reordering, "amd", "Determines what sort of preordering will be used"
								 " on the matrix. Choices are 'amd', 'rcm', and 'none'.");
								 
DEFINE_string(equil, "bunch", "Decides if the matrix should be equilibriated (in the max-norm) "
						 "before factoring is done. Choices are 'none', 'bunch', and 'iter'");
						 
DEFINE_bool(save, true, "If yes, saves the factors (in matrix-market format) into a folder"
						 "called output_matrices/ in the same directory as ldl_driver.");

DEFINE_bool(display, false, "If yes, outputs a human readable version of the factors onto"
							" standard out. Generates a large amount of output if the "
							"matrix is big.");
		

int main(int argc, char* argv[])
{
	std::string usage("Performs an incomplete LDL factorization of a given matrix.\n"
				      "Sample usage:\n"
					  "\t./skew_ldl_driver -filename=test_matrices/skew_testmat1.mtx "
									 "-fill=2.0 -display=true -save=false\n"
					  "Additionally, these flags can be loaded from a single file "
					  "with the option -flagfile=[filename].");
				 
	google::SetUsageMessage(usage);
	google::ParseCommandLineFlags(&argc, &argv, true);
	
	if (FLAGS_filename.empty()) {
		std::cerr << "No file specified! Type ./skew_ldl_driver --help for a description of the program parameters." << std::endl;
		return 0;
	}
	
	skew_solver<double> solv;
	solv.load(FLAGS_filename);
	
	//default reordering scheme is AMD
	solv.set_reorder_scheme(FLAGS_reordering.c_str());
	
	//default is equil on
	solv.set_equil(FLAGS_equil.c_str()); 
	
	solv.solve_pp(FLAGS_fill, FLAGS_tol, FLAGS_pp_tol);
	
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

