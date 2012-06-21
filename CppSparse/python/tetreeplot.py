import cppsparse as c
trips = [
    (1,6),(1,7),
    (2,3),(2,8),
    (3,10),(3,11),
    (4,6),(4,10),
    (5,8),(5,11),
    (6,9),(6,10),
    (7,11),
    (8,10),(8,11),
    (10,11)]

def build_matrix():
    N = 12;
    t = c.dtrp(N, N ,4*N)
    for row in range(N):
        t.push_back(int(row), int(row), 5.0);

    for row, col in trips:
        t.push_back(int(row), int(col), 2.0);
        t.push_back(int(col), int(row), 2.0);

    return c.dcsc(t)

if __name__=="__main__":
    s = build_matrix();
    s.etree_plot("elim_tree.dot")

    
