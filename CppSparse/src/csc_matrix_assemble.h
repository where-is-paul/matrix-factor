// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_ASSEMBLE_H_
#define _CSC_MATRIX_ASSEMBLE_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: assemble ()
{
    static const idx_type invalid_idx = static_cast<idx_type> (-1);
    idx_vector_type work(m_n_rows, invalid_idx);

    idx_type current_offset = 0;

    for (idx_type j = 0; j < m_n_cols; j ++)
    {
        /**
         * Save the start index of the current column. We update
         * m_col_idx[j] _after_ we have processed the j-th column.
         */

        idx_type current_col_start = current_offset;

        for (idx_type offset = m_col_idx[j]; offset < m_col_idx[j + 1]; offset ++)
        {
            idx_type row = m_row_idx[offset];
            el_type  val = m_x [offset];
            if (val > 0)
            {
                /**
                 *  This is the standard "trick" of using a work array
                 *  to compute the row offsets.  Initially, all
                 *  elements in the work array are filled with an
                 *  invalid index. Then, for each row in the matrix,
                 *  we see if the work element is either invalid or
                 *  belongs to a previous column. If the row has
                 *  already occured in the current column, then we
                 *  simply sum it to the existing value at the
                 *  offset. Otherwise, we update the work with the
                 *  current offset.
                 */
                if (work[row] == invalid_idx || work[row] < current_col_start)
                {
                    work[row] = current_offset;
                    m_row_idx[current_offset] = row;
                    m_x      [current_offset] = val;
                    current_offset ++;
                }
                else
                {
                    m_x [work[row]] += val;
                }
            }
        }

        /**
         *  At this point, we have updated the j-th column, so it is
         *  safe to update the start of the column based on the saved
         *  column start.
         */
        m_col_idx[j] = current_col_start;
    }
    m_col_idx.back() = current_offset;
    m_row_idx.resize(current_offset);
    m_x.resize      (current_offset);
    
    /**
     * Now, we need to sort the row indices
     */
    sort();
}

#endif
