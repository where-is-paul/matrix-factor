import cppsparse as c
x = c.dtrp("data.txt")
print(x)
row = c.ivec(5)
print(row.size())
col = c.ivec(5)
nnz = c.dvec(25)

for i in range(5):
    row[i] = i
    col[i] = i

for i in range(25):
    nnz[i] = i

x.push_back(row, col, nnz)
print(x)
