// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_TRIANGULAR_H_
#define _CSC_MATRIX_TRIANGULAR_H_

template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: triangular (uplo_type uplo) const
{
    if (uplo == UPPER_TRIANGULAR)
    {
        return filter([](idx_type row, idx_type col, el_type val) {return col >= row;});
    }
    else
    {
        return filter([](idx_type row, idx_type col, el_type val) {return col <= row;});
    }
}
#endif
