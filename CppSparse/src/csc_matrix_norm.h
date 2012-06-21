// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_NORM_H_
#define _CSC_MATRIX_NORM_H_

template<class idx_type, class el_type>
el_type csc_matrix<idx_type, el_type> :: norm_one () const
{
    el_type max_val = 0.0;

    for (idx_type j = 0; j < m_n_cols; j ++)
    {
        el_type sum = std::accumulate(m_x.begin() + m_col_idx[j], m_x.begin() + m_col_idx[j+1], 0.0, 
                                      [](el_type a, el_type b)  { return std::abs(a) + std::abs(b); });
        max_val = std::max(sum, max_val);
    }

    return max_val;
}


template<class idx_type, class el_type>
el_type csc_matrix<idx_type, el_type> :: norm_inf () const
{

    elt_vector_type sum (n_rows());

    for (idx_type el = 0; el < nnz(); el ++)
    {
        sum[m_row_idx[el]] += std::abs(m_x[el]);
    }
    
    return std::accumulate(sum.begin(), sum.end(), 0.0, [](el_type a, el_type b){return std::max(a, b);});
}

template<class idx_type, class el_type>
el_type csc_matrix<idx_type, el_type> :: norm_fro () const
{
    el_type sum = std::accumulate(m_x.begin(), m_x.end(),  0.0, 
                                  [](el_type a, el_type b) { return a*a + b*b; });
    return sqrtf(sum);
}


template<class idx_type, class el_type>
el_type csc_matrix<idx_type, el_type> :: norm (norm_type nt) const
{
    switch (nt)
    {
    case NORM_ONE:
        return norm_one();
    case NORM_INF:
        return norm_inf();
    case NORM_FRO:
        return norm_fro();
    default:
        throw std::logic_error("Norm must be ONE, INF, or FRO");
    }
}

#endif
