import cppsparse
import numpyinterface as ni
import numpy
import unittest

class TestDenseTriangularSolve(unittest.TestCase):

    def setUp(self):
        self.data = [(10, 0.5), (20, 0.4), (40, 0.5), (80, 0.4)]

    def verifySolve(self, A, b, trans, xAct):
        m = A.n_rows()
        res = cppsparse.dvec(m, 0.0);

        if trans:
            normA = A.norm(cppsparse.NORM_INF)
            A.gatxpy(xAct, res)
        else:
            normA = A.norm(cppsparse.NORM_ONE)
            A.gaxpy(xAct, res)

        normRes = numpy.linalg.norm(ni.nvec(b) - ni.nvec(res), numpy.inf)
        self.assertTrue(normRes/(numpy.finfo('double').eps * normA * m) < 100, 
                        'Residual norm = %f'%(normRes/(numpy.finfo('double').eps * normA * m)))
        
    def testLSolve(self):   
        for m, s in self.data:
            L = ni.lower_tri_sp(m, m, s)
            b = ni.vect(ni.np.arange(m, dtype=ni.np.float64))
            x = cppsparse.dvec(b)
            L.lsolve(x)
            self.verifySolve(L, b, False, x)

    
    def testLTSolve(self):
        for m, s in self.data:
            L = ni.lower_tri_sp(m, m, s)
            b = ni.vect(ni.np.arange(m, dtype=ni.np.float64))
            x = cppsparse.dvec(b)
            L.ltsolve(x)
            self.verifySolve(L, b, True, x)


    def testUSolve(self):
        for m, s in self.data:
            U = ni.upper_tri_sp(m, m, s)
            b = ni.vect(ni.np.arange(m, dtype=ni.np.float64))
            x = cppsparse.dvec(b)
            U.usolve(x)
            self.verifySolve(U, b, False, x)

    def testUTSolve(self):
        for m, s in self.data:
            U = ni.upper_tri_sp(m, m, s)
            b = ni.vect(ni.np.arange(m, dtype=ni.np.float64))
            x = cppsparse.dvec(b)
            U.utsolve(x)
            self.verifySolve(U, b, True, x)

    
if __name__=="__main__":
    unittest.main()

