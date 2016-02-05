include ../make_defs.inc
ILDL_TARGET = ildl

MEXFLAGS = -largeArrayDims 
MEXEXTRA = utils

SRCDIR := ..
INC_MC64 = $(addprefix -I, $(SRCDIR)/include/hsl_mc64)
LIB_MC64 = $(SRCDIR)/lib/libhsl_mc64.a
LIB_GFORTRAN = $(SRCDIR)/lib/libhsl_mc64.a

all: $(addsuffix .$(MEX_EXT), $(ILDL_TARGET))

$(addsuffix .$(MEX_EXT), $(ILDL_TARGET)): $(addsuffix .cpp, $(ILDL_TARGET) $(MEXEXTRA)) $(LIB_MC64)
	$(MATLAB_BIN) $(INC_MC64) $(MEXFLAGS) $(LINKFLAGS) $?
    
clean:
	$(RM) $(addsuffix .$(MEX_EXT), $(ILDL_TARGET))
    
test: 
	matlab -nodesktop -nosplash -r 'test_mex; exit'

.PHONY : all clean test