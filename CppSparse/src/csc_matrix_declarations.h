// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_DECLARATIONS_H_
#define _CSC_MATRIX_DECLARATIONS_H_

#include "abstract_sparse_matrix.h"

template<class idx_type, class el_type> 
class csc_matrix : public abstract_sparse_matrix<idx_type, el_type>
{
public:
    typedef csc_matrix<idx_type, el_type> csc_t;

    using abstract_sparse_matrix<idx_type, el_type>::m_row_idx;
    using abstract_sparse_matrix<idx_type, el_type>::m_col_idx;
    using abstract_sparse_matrix<idx_type, el_type>::m_x;
    using abstract_sparse_matrix<idx_type, el_type>::m_n_rows;
    using abstract_sparse_matrix<idx_type, el_type>::m_n_cols;
    using abstract_sparse_matrix<idx_type, el_type>::n_rows;
    using abstract_sparse_matrix<idx_type, el_type>::n_cols;
    using abstract_sparse_matrix<idx_type, el_type>::shape;


    typedef typename abstract_sparse_matrix<idx_type, el_type>::idx_vector_type idx_vector_type;
    typedef typename abstract_sparse_matrix<idx_type, el_type>::elt_vector_type elt_vector_type;

    
private:
    static const idx_type npos;

    void compress(const triplet_matrix<idx_type, el_type>& triplet);    
    
    template<class predicate>
    void keep (const predicate& pred);

    template<class predicate>
    csc_t filter (const predicate& pred) const;

    /// Gathers the contents of a dense vector based on 
    /// non-zero pattern of the matrix for a given column
    void gather (idx_type column, const elt_vector_type& vec);
    
    /// Compute x += A(:, k) * beta and updates the non-zero position of
    /// the result based on the non-zero values of x
    void scatter (csc_t& result, idx_type column, el_type beta, 
                  idx_vector_type& work, elt_vector_type& x, idx_type mark, idx_type& nnz) const;

    el_type norm_one () const;
    el_type norm_inf () const;
    el_type norm_fro () const;

    csc_t col_permute (const idx_vector_type& vec) const;
    csc_t row_permute (const idx_vector_type& vec) const;
    csc_t sym_permute (const idx_vector_type& vec) const
    {
        throw std::runtime_error("Not implemented.");
    }

    el_type dot_is_sorted (const csc_matrix<idx_type, el_type>& y) const;
    el_type dot_no_sorted (const csc_matrix<idx_type, el_type>& y) const;

    void do_hcat(idx_type col_offset, const csc_matrix<idx_type, el_type>& x);
    void do_vcat(idx_type row_offset, std::vector<idx_type>& work, const csc_matrix<idx_type, el_type>& x);

    friend csc_matrix<idx_type, el_type> hcat (const csc_matrix<idx_type, el_type>& x, 
                                               const csc_matrix<idx_type, el_type>& y)
    {
        if (x.n_rows() != y.n_rows())
        {
            throw std::logic_error("Number of rows is inconsistent");
        }
        
        idx_type result_nrows = x.n_rows();
        idx_type result_ncols = x.n_cols() + y.n_cols();
        idx_type result_nnz   = x.nnz() + y.nnz();
        csc_matrix<idx_type, el_type> result(result_nrows, result_ncols, result_nnz);
        result.m_row_idx.resize(result_nnz);
        result.m_x.resize(result_nnz);
        result.do_hcat(0, x);
        result.do_hcat(x.n_cols(), y);
        return result;
    }

    friend csc_matrix<idx_type, el_type> vcat (const csc_matrix<idx_type, el_type>& x,
                                               const csc_matrix<idx_type, el_type>& y)
    {
        if (x.n_cols() != y.n_cols())
        {
            throw std::logic_error("Number of columns is inconsistent");
        }
        
        idx_type result_nrows = x.n_rows() + y.n_rows();
        idx_type result_ncols = x.n_cols();
        idx_type result_nnz   = x.nnz() + y.nnz();
        csc_matrix<idx_type, el_type> result(result_nrows, result_ncols, result_nnz);

        // The column indices for the result is simply the sum of
        // the column indices of the individual arrays to be
        // concatenated
        std::transform(x.m_col_idx.begin(), x.m_col_idx.end(), y.m_col_idx.begin(), result.m_col_idx.begin(),
                       std::plus<idx_type>());
            
        // VCAT requires a work vector of the same size as the
        // number of columns of X and Y. It initially stores the
        // column indices, just like "csc_matrix_compress"
        std::vector<idx_type> work(result.m_n_cols);
        std::copy (result.m_col_idx.begin(), result.m_col_idx.end()-1, work.begin());

        result.m_row_idx.resize(result_nnz);
        result.m_x.resize(result_nnz);
        result.do_vcat(0, work, x);
        result.do_vcat(x.n_rows(), work, y);
        return result;           
    }

    /**
     *  For CSC matrix, squeeze out zero values, sum up duplicate
     *  elements and sort the row indices.
     */
    void assemble ();

    void dfs(idx_type col, idx_vector_type& xi, 
             idx_vector_type& pstack, const idx_vector_type& pinv, idx_type& top) const;


    /**
     * Symbolic analysis of a permuted lower-triangular matrix.
     * Exercise 3.3
     */
    void lower_row_perm (idx_vector_type& perm) const;    

    /**
     * Symbolic analysis of a permuted upper-triangular matrix.
     * Exercise 3.4
     */
    void upper_row_perm (idx_vector_type& perm) const;

    void lower_col_perm (idx_vector_type& perm) const;

    void upper_col_perm (idx_vector_type& perm) const;

    /**
     * Method to construct a postordering for a given column.
     * @param col The parent whos children must be enumerated.
     * @param first_child The list of all children. This is mutated inside this method.
     * @param next_sibling The list of all first siblings.
     * @param post The post-ordered tree, updated with the information for col.
     * @param stack The recursion stack, must be emptied out for each column.
     */ 
    void dfs_tree (idx_type col, idx_vector_type& first_child, 
                   const idx_vector_type& next_sibling, 
                   idx_vector_type& post, 
                   idx_vector_type& stack) const;

public:

    csc_matrix (idx_type n_rows = 0, idx_type n_cols = 0, idx_type nz_max = 0): 
        abstract_sparse_matrix<idx_type, el_type> (n_rows, n_cols, nz_max) 
    {
        m_col_idx.resize (n_cols + 1);
    }


    csc_matrix (const triplet_matrix<idx_type, el_type>& triplet, bool assemble = true) : 
        abstract_sparse_matrix<idx_type, el_type> (triplet.n_rows(), triplet.n_cols())
    {
        compress(triplet);
        if (assemble)
        {
            this -> assemble();
        }
    }

    /*
      csc_matrix (csc_matrix&& other)
      {
      *this = std::move (other);
      }

      csc_matrix& operator= (csc_matrix&& other) 
      {    
      if (this != &other)
      {
      m_row_idx  = other.row_idx;
      m_col_idx  = other.col_idx;
      m_x        = other.m_x;
      }
      return *this;
      }
    */

    /// In compressed sparse column storage, the col_idx array is of size N + 1
    /// 
    /// col_idx[j] gives the starting position of the first non-zero element in column j
    /// 
    /// Hence col_idx[j+1] - col_idx[j] gives the total number of
    /// non-zero values in column j and therefore, col_idx[N] gives the
    /// total number of non-zero elements in the matrix.
    /// 
    /// row_idx[j] and nz_val[j] are arrays of size NNZ, so col_idx[N] == row_idx.size()  
    virtual idx_type nnz() const
    {
        return m_row_idx.size();
    }  

    /**
     * Returns a string representation of a CSC matrix.
     * @returns A string representation of this matrix.
     */
    std::string to_string () const;

    /**
     * Returns a graph representation of this matrix
     * as a DOT-formatted directed graph. 
     * @stream The output stream to which the graph must be written to.
     * @trans If True, output the graph of A<sup>T</sup>.
     */
    void        to_graph (std::ostream& stream, bool trans = false) const;
    void        to_graph (const std::string& file_name, bool trans = false) const;
    std::string to_graph (bool trans = false) const;

    /**
     * Sort the row indices in each column, c.f. Exercises 2.7 and 2.11.
     */
    void sort ();

    /**
     * Returns a triplet representation of this matrix.
     * @return A triplet matrix consisting of uncompressed column indices.
     */
    triplet_matrix<idx_type, el_type> find () const;

    /**
     * Computes y = A * x + y
     * @param x The vector <em>x</em>
     * @param y The vector <em>y</em>
     */
    void gaxpy (const el_type* x, el_type* y) const;

    void gaxpy (const std::vector<el_type>& x, std::vector<el_type>& y) const;

    /**
     * Computes y = A<sup>T</sup> * x + y.
     * @param x The vector <em>x</em>.
     * @param y The vector <em>y</em>.
     */
    void gatxpy (const std::vector<el_type>& x, std::vector<el_type>& y) const;

    /**
     * Computes y = A * x + y when A is symmetric. Only one half of A is accessed.
     * @param x The vector <em>x</em>.
     * @param y The vector <em>y</em>.
     * @param side If side = UPPER_TRIANGULAR, only the
     * upper-triangular part of A is accessed otherwise the
     * lower-triangular part is accessed.
     */
    void sym_gaxpy (const el_type* x, el_type* y, uplo_type side = UPPER_TRIANGULAR) const;

    void sym_gaxpy (const std::vector<el_type>& x, std::vector<el_type>& y, uplo_type side = UPPER_TRIANGULAR) const;

    /**
     * Computes the dot product of two vectors, x.y
     * @param y The other sparse vector.
     * @param is_sorted If true, assume columns of x and y are sorted.
     */
    el_type dot (const csc_matrix<idx_type, el_type>& y, bool is_sorted = true) const;

    /**
     * Compute RxAxC where R and C are diagonal matrices.
     * @param r The diagonal matrix to premultiply.
     * @param c The diagonal matrix to postmultiply.
     */
    void scale (const std::vector<el_type>& r, const std::vector<el_type>& c);

    void keep ()
    {
        keep ([](idx_type i, idx_type j, el_type v) {return v != 0;});
    }

    void transpose(csc_matrix<idx_type, el_type>& result) const;

    csc_t transpose() const;

    csc_t multiply (const csc_t& b) const;

    csc_t multiply2 (const csc_t& b) const;

    csc_t add (const csc_t& b, el_type alpha = 1.0, el_type beta = 1.0) const;

    el_type norm (norm_type type = NORM_ONE) const;
    
    csc_t permute (const idx_vector_type& p, perm_type type = PERM_COL) const;

    csc_t permute (const idx_vector_type& p, const idx_vector_type& q) const;   

    csc_t triangular (uplo_type = UPPER_TRIANGULAR) const;

    virtual std::vector<el_type> full () const;

    void sum_duplicates ();

    /**
     * Comparison operator for CSC matrices.
     * @param other The matrix to be compared
     * @return True if the two matrices are the same under the following assumptions:
     * - They have the same number of non-zero entries (i.e., all numerical zero entries have been removed)
     * - Each column consists of the same row indices in sorted order
     * - The non-zero values in column are the same
     */
    bool operator==(const csc_matrix<idx_type, el_type>& other) const;

    bool is_symmetrix () const
    {
        return (*this) == transpose();
    }
    
    csc_t subsref_range(idx_type i1, idx_type i2, idx_type j1, idx_type j2) const;

    csc_t subsref_vector (const idx_vector_type& rows, const idx_vector_type& cols) const;

    void lsolve  (elt_vector_type& b) const;

    void ltsolve (elt_vector_type& b) const;

    void usolve  (elt_vector_type& b) const;

    void utsolve (elt_vector_type& b) const;

    void lsolve (const csc_t& b, idx_vector_type& x, idx_vector_type& xi, idx_type k = 0) const;

    void lsolve(const csc_t& b, idx_vector_type& x, idx_type k = 0) const
    {
        idx_vector_type xi(m_n_cols);
        lsolve(b, x, xi, k);
    }

    idx_type reach(const csc_t& b, idx_type k, const idx_vector_type& pinv, idx_vector_type& xi) const;

    void lower_row_perm_solve (elt_vector_type& b) const;

    void upper_row_perm_solve (elt_vector_type& b) const;

    void lower_col_perm_solve (elt_vector_type& b) const;

    void upper_col_perm_solve (elt_vector_type& b) const;
    
    /**
     * Solve Ax = b when A is a permuted lower-triangular matrix and b is a dense vector.
     * @param b The right-hand side to solve against. It is overwritten with the solution.
     */
    void triangular_solve (elt_vector_type& b) const;

    idx_vector_type etree(bool trans=false) const;

    void        etree_plot (std::ostream& stream, bool trans = false) const;
    void        etree_plot (const std::string& file_name, bool trans = false) const;
    std::string etree_plot (bool trans = false) const;

    idx_type ereach (const idx_vector_type& parent, 
                     idx_type col, idx_vector_type& stack, idx_vector_type& work) const;

    idx_vector_type ereach(idx_type col) const;

    /**
     * Postorder the elimination tree for a symmetrix matrix.
     * @param parent The array representing the elimination tree for this matrix.
     * @returns The postordered tree such that the <em>d</em>
     * decendents of a verted <em>k</em> are numbered consecutively.
     */
    idx_vector_type postorder(const idx_vector_type& parent) const;
};

template<class idx_type, class el_type>
const idx_type csc_matrix<idx_type, el_type>::npos = static_cast<idx_type>(-1);

#endif
