include ../make_defs.inc
MEX_EXT = $(shell $(MATLAB_BIN)/mexext)

ILDL_TARGET = ildl

MEXFLAGS = -largeArrayDims 
MEXEXTRA = utils

all: $(addsuffix .$(MEX_EXT), $(ILDL_TARGET))

$(addsuffix .$(MEX_EXT), $(ILDL_TARGET)): $(addsuffix .cpp, $(ILDL_TARGET) $(MEXEXTRA))
	$(MATLAB_BIN)/mex $(MEXFLAGS) $?
    
clean:
	$(RM) $(addsuffix .$(MEX_EXT), $(ILDL_TARGET))
    
test: 
	matlab -nodesktop -nosplash -r 'test_mex; exit'

.PHONY : all clean test
