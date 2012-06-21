import cppsparse as c
import numpy
    
def nfull(x):
    m, n = (x.n_rows(), x.n_cols())
    y = numpy.array(x.full(), order = 'F')
    y.shape = m, n
    return y.transpose()

def nvec(x):
    return numpy.array(x)

def random_col_vector(m, nnz):
    irn = numpy.random.random_integers(0, m-1, nnz).astype(numpy.uint64)
    irn = numpy.random.permutation(numpy.unique(irn));
    nzv = numpy.random.rand(nnz)    
    triplet_matrix = c.dtrp(m, 1, nnz)

    for ir, nz in zip(irn, nzv):
        triplet_matrix.push_back(int(ir), 0, nz)

    return c.dcsc(triplet_matrix)

x = random_col_vector(2000, 100)
y = random_col_vector(2000, 100)

print x.dot(y, False)
print numpy.dot(nfull(x), nfull(y).T)[0][0]
