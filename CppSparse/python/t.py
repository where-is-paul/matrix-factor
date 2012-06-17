import cppsparse as c
x = c.dtrp("data.txt")
print(x)
y = c.dcsc(x)
print(y)
print(y == y)
t = y.transpose()
print(t)
u = t.transpose()
print(u)
v = u.transpose()
print(v.find())


