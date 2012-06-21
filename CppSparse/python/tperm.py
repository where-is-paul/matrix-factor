import numpyinterface as n
a = n.c.dtrp(3,3)
i = [0,1,2,0,2,0]
j = [0,0,0,1,1,2]
v = [1,2,3,4,5,6]
map(lambda i, j, v: a.push_back(int(i), int(j), float(v)), i, j, v)
a = n.c.dcsc(a)
print(n.full(a))
b = n.c.dvec(3)
a.triangular_solve(b)

