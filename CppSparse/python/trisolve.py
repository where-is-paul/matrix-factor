import numpyinterface as n
c=n.c
t = c.dtrp()
t.push_back(2,0,1.0);
t.push_back(0,1,2.0);
t.push_back(1,1,3.0);
t.push_back(2,1,4.0);
t.push_back(1,2,5.0);
t.push_back(2,2,6.0);
x = c.dcsc(t)
b = c.dvec(3)
x.triangular_solve(b)

t = c.dtrp()
t.push_back(0,0,1.0)
t.push_back(1,0,2.0)
t.push_back(0,1,3.0)
t.push_back(1,1,4.0)
t.push_back(2,1,5.0)
t.push_back(0,2,6.0)
x = c.dcsc(t)
x.triangular_solve(b)

