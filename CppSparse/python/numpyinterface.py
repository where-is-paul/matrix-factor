import cppsparse as c
import numpy as np

def full(x):
    m, n = (x.n_rows(), x.n_cols())
    y = np.array(x.full(), order = 'F')
    if m > 1 and n > 1:
        y.shape = n, m
        return y.transpose()
    else:
        y.shape = m, n
        return y

def nvec(x):
    return np.array(x)

def spcolvec(x):
    y = c.dtrp(len(x), 1)
    for i in range(len(x)):
        y.push_back(i, 0, float(x[i]))
    z = c.dcsc(y)
    return z

def vect(n):
    v = None
    if n.dtype == np.float64:
        v = c.dvec();
        for elem in n:
            v.push_back(float(elem))
    elif n.dtype == np.int32 or n.dtype == np.int64:        
        v = c.ivec()    
        for elem in n:
            v.push_back(int(elem))
    else:
        raise RuntimeError("Array must be float or integer")
    return v

def lower_tri_full(m, n):
    l = c.dtrp(m, n)
    for col in range(n):
        for row in range(col, m):
            l.push_back(row, col, float(row + col * m))
    return c.dcsc(l)

def upper_tri_full(m, n):
    u = c.dtrp(m, n)
    for col in range(n):
        for row in range(col, m):
            u.push_back(col, row, float(row + col * m))
    return c.dcsc(u)

def sprand(m, n, s):
    nnz = int(m * n * s)
    row_idx = vect(np.random.random_integers(0, m - 1, nnz))
    col_idx = vect(np.random.random_integers(0, n - 1, nnz))
    vals    = vect(np.random.rand(nnz))
    
    for d in range(0, min(m, n)):
        row_idx.push_back(int(d))
        col_idx.push_back(int(d))
        vals.push_back(float(1.0))

    tr = c.dtrp(m, n)

    for i, j, v in zip(row_idx, col_idx, vals):
        tr.push_back(i, j, v)
    
    return c.dcsc(tr)
    
def upper_tri_sp(m, n, s = 0.2):
    return sprand(m, n, s).triangular(c.UPPER_TRIANGULAR)


def lower_tri_sp(m, n, s = 0.2):
    return sprand(m, n, s).triangular(c.LOWER_TRIANGULAR)
