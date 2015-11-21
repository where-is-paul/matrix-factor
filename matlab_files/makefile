CXXFLAGS = -Wall -O3 -std=c++0x

FILES = utils.cpp ../*.cpp ../*.h

all : ildl.mex*

UNAME := $(shell uname)
MATLABDIR := /Applications/Matlab/MATLAB_R2015b.app

ifeq ($(UNAME), Linux)
MEXNAME = mex
else ifeq ($(UNAME), Darwin)
MEXNAME = mex
else
MEXNAME = mex.bat
endif

MEX = $(MATLABDIR)/bin/$(MEXNAME)

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
