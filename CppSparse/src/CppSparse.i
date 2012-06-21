%module cppsparse
%include "std_string.i"
%include "std_vector.i"

%{
// Includes the header in the wrapper code 
#include "triplet_matrix.h"
#include "csc_matrix.h"
#include <sstream>
#include <iostream>
%}
 

/* Parse the header file to generate wrappers */
%include "abstract_sparse_matrix.h"
%include "triplet_matrix.h"
%include "csc_matrix.h"
%include "csc_matrix_declarations.h"

%template(dabs) abstract_sparse_matrix<size_t, double>;
%template(dcsc) csc_matrix<size_t, double>;
%template(dtrp) triplet_matrix<size_t, double>;
%template(dvec) std::vector<double>;
%template(ivec) std::vector<size_t>;

