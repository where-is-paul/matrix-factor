CXXFLAGS = -Wall -O3 -std=c++0x

FILES = utils.cpp ../*.cpp ../*.h

all : ildl.mex*

UNAME := $(shell uname)

# The part below requires some set up. If this makefile doesn't work out of the box,
# you'll need to find the mex script in your MATLAB folder below
ifeq ($(UNAME), Linux)
    MEX = mex
else ifeq ($(UNAME), Darwin)
    MEXNAME = mex
    MATLABDIR := /Applications/Matlab/MATLAB_R2015b.app
    MEX = $(MATLABDIR)/bin/$(MEXNAME)
else
    MEXNAME = mex.bat
    MATLABDIR := C:/Program\ Files/MATLAB/MATLAB\ Production\ Server/R2015a
    MEX = $(MATLABDIR)/bin/$(MEXNAME)
endif

MEXFLAGS = -largeArrayDims
MEXEXTRA = ildl.cpp utils.cpp

%.mex*: %.cpp
	$(MEX) $(MEXFLAGS) $(MEXEXTRA)

ildl.mex* : ildl.cpp
    
clean:
	$(RM) *.mex*
    
test: 
	matlab -nodesktop -nosplash -r test_mex

debug: 
	matlab -Dgdb
.PHONY : all clean test
