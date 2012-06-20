import cppsparse as c
x = c.dtrp("data.txt")
print(x)
y = c.dcsc(x)
z = c.vcat(y, y)
print(z)
print(z.find())
