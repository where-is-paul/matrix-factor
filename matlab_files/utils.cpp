#include <algorithm>
#include "utils.h"

using namespace std;

bool is_double_scalar(const mxArray* a) {
    return mxIsDouble(a)
            && !mxIsComplex(a)
            && mxGetNumberOfElements(a) == 1;
}

bool is_vector(const mxArray* a) {
    const mwSize* dims = mxGetDimensions(a);
    return mxGetNumberOfDimensions(a) == 2
            && min(dims[0], dims[1]) <= 1;
}

bool is_double_vector(const mxArray* a) {
    return mxIsDouble(a)
            && !mxIsComplex(a)
            && is_vector(a);
}

bool is_string(const mxArray* a) {
    return mxIsChar(a) && is_vector(a);
}

double parse_double(const mxArray* raw_double) {
    if (!is_double_scalar(raw_double))
        mexErrMsgTxt("Argument must be a real scalar.");
    double res = *((double*) mxGetData(raw_double));

    return res;
}