import cppsparse as c

def nfull(x):
    import numpy 
    m, n = (x.n_rows(), x.n_cols())
    y = numpy.array(x.full(), order = 'F')
    y.shape = m, n
    return y.transpose()

def nvec(x):
    import numpy
    return numpy.array(x)

lower_tri = c.dtrp(5, 5)
full_sym  = c.dtrp(5, 5)
upper_tri = c.dtrp(5, 5)

for col in range(5):
    for row in range(col, 5):
        lower_tri.push_back(row, col, float(row + col * 5.0))
        full_sym.push_back (row, col, float(row + col * 5.0))
        upper_tri.push_back (col, row, float(row + col * 5.0))
        if row != col:
            full_sym.push_back (col, row, float(row + col * 5.0))


lower_tri = c.dcsc(lower_tri)
upper_tri = c.dcsc(upper_tri)
full_sym  = c.dcsc(full_sym)

x = c.dvec(range(5))

y = c.dvec(5)
lower_tri.sym_gaxpy(x, y, c.LOWER_TRIANGULAR)
print(nvec(y))

y = c.dvec(5)
upper_tri.sym_gaxpy(x, y, c.UPPER_TRIANGULAR)
print(nvec(y))

y = c.dvec(5)
full_sym.gaxpy(x, y)
print(nvec(y))



