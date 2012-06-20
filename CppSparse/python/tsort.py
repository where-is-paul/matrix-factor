import cppsparse as c

def random_triplet (m, n, nnz):
    import numpy as np
    row_idx = np.random.random_integers(0, m-1, nnz)
    col_idx = np.random.random_integers(0, n-1, nnz)
    nz__val = np.random.rand(10)

    trip = c.dtrp(m, n, nnz)
    
    for i,j,v in zip(row_idx, col_idx, nz__val):
        print i, j, v
        trip.push_back(int(i), int(j), float(v))

    return trip


x = c.dcsc(random_triplet(8, 8, 20))
##print(x)
#print(x.find())
x.sort()
##print(x.find())








