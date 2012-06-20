import cppsparse as cp
x = cp.dcsc(cp.dtrp("fem.dat"))
x.sum_duplicates()
print(x)
