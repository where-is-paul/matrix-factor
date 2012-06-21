def patchfile(fname):
    fnamepatched = fname + '.p'
    f = open(fname)
    g = open(fnamepatched,'w')

    for line in f:
        if '#include <Python.h>' in line:
            g.write('#include <cmath>\n')
            g.write(line)
        else:
            g.write(line)

    g.flush()
    f.close()
    g.close()
    import os
    os.unlink(fname)
    os.rename(fnamepatched, fname)
    

if __name__== "__main__":
    import sys
    patchfile(sys.argv[1])
        
