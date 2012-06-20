// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_H_
#define _CSC_MATRIX_H_

#include "abstract_sparse_matrix.h"

template<class idx_type, class el_type>
class triplet_matrix;

enum norm_type {NORM_ONE, NORM_INF, NORM_FRO};
enum perm_type {PERM_ROW, PERM_COL, PERM_SYM};
enum uplo_type {UPPER_TRIANGULAR, LOWER_TRIANGULAR};
template<class idx_type, class el_type> 
class csc_matrix;

#include "csc_matrix_declarations.h"
#include "csc_matrix_compress.h"
#include "csc_matrix_assemble.h"
#include "csc_matrix_find.h"
#include "csc_matrix_sort.h"
#include "csc_matrix_scale.h"
#include "csc_matrix_to_string.h"
#include "csc_matrix_to_graph.h"
#include "csc_matrix_gaxpy.h"
#include "csc_matrix_gatxpy.h"
#include "csc_matrix_sym_gaxpy.h"
#include "csc_matrix_dot.h"
#include "csc_matrix_compress.h"
#include "csc_matrix_transpose.h"
#include "csc_matrix_sum_duplicates.h"
#include "csc_matrix_keep.h"
#include "csc_matrix_scatter.h"
#include "csc_matrix_gather.h"
#include "csc_matrix_multiply.h"
#include "csc_matrix_multiply2.h"
#include "csc_matrix_add.h"
#include "csc_matrix_norm.h"
#include "csc_matrix_perm.h"
#include "csc_matrix_full.h"
#include "csc_matrix_operator_eq.h"
#include "csc_matrix_hcat.h"
#include "csc_matrix_vcat.h"
#include "csc_matrix_subsref.h"
#include "csc_matrix_filter.h"
#include "csc_matrix_triangular.h"
#include "csc_matrix_triangular_solve.h"
#include "csc_matrix_reach.h"
#include "csc_matrix_sp_triangular_solve.h"
#include "csc_matrix_perm_row_triangular_solve.h"
#include "csc_matrix_perm_col_triangular_solve.h"
#include "csc_matrix_perm_triangular_solve_general.h"
#include "csc_matrix_etree.h"
#include "csc_matrix_ereach.h"
#include "csc_matrix_postorder.h"
#endif 
