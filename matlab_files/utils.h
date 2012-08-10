#ifndef _MEX_UTILS
#define _MEX_UTILS
#include "../source/lilc_matrix.h"
#include <mex.h>

// validation methods
bool is_double_scalar(const mxArray* a);
bool is_vector(const mxArray* a);
bool is_double_vector(const mxArray* a);
bool is_string(const mxArray* a);

// parsing methods
double parse_double(const mxArray* raw_alpha);

#endif