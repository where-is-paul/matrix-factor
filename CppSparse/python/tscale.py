import cppsparse as c
x = c.dtrp("fem.dat")
y = c.dcsc(x);
y.sum_duplicates()
m = y.n_rows()
n = y.n_cols()
print m, n
cs = c.dvec(n, 2.0)
rs = c.dvec(m, 3.0)
print(y.find())
y.scale(rs, cs)
print(y.find())



