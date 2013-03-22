CXXFLAGS = -Wall -O3 -std=c++0x

FILES = utils.cpp ../*.cpp ../*.h

all : ildl.mex*

UNAME := $(shell uname)
MEX = mex

ifeq ($(UNAME), Linux)
MEX = mex
else
MEX = mex.bat
endif

MEXFLAGS =-largeArrayDims
MEXEXTRA =ildl.cpp utils.cpp

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