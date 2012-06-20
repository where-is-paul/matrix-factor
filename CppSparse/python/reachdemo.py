import cppsparse as c
##import numpyinterface as n

m = 15
conn = [(1,3), (3,5), (3,8), (5,8), (8,10), (6,9), (6,10), (2,4), (2,7), (2,9), (4,9),
        (6,9), (6,10), (7,10), (7,11), (8,10), (9,12), (9,13), (10,11), (10,13), (10,14), (11,12), (11,13), (12,13), (13,14)]

trip = c.dtrp(m, m, 2*m)

for i in range(m):
    trip.push_back(i, i, 1.0)

for i, j in conn:
    trip.push_back(i, j, 1.0)

G = c.dcsc(trip).transpose();
print(G)
#G.to_graph("tim.dot", True)

rhs = c.dtrp(15,1,4)
rhs.push_back(4,0,1.0)
rhs.push_back(6,0,1.0)

B = c.dcsc(rhs)
top = G.reach(B, 0, c.ivec(m), c.ivec(m));


