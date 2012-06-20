import cppsparse as c
x = c.dtrp("fem.dat")
print(x)
y = c.dcsc(x);
y.sum_duplicates()
print(y)
print(y.find())
