include ../make_defs.inc
ILDL_TARGET = ildl

MEXFLAGS = -largeArrayDims 
MEXEXTRA = utils

SRCDIR := ..

all: $(addsuffix .$(MEX_EXT), $(ILDL_TARGET))

$(addsuffix .$(MEX_EXT), $(ILDL_TARGET)): $(addsuffix .cpp, $(ILDL_TARGET) $(MEXEXTRA))
	$(MATLAB_BIN) $(MEXFLAGS) $(LINKFLAGS) $?
    
clean:
	$(RM) $(addsuffix .$(MEX_EXT), $(ILDL_TARGET))
    
test: 
	matlab -nodesktop -nosplash -r 'test_mex; exit'

.PHONY : all clean test