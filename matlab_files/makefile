CXXFLAGS = -Wall -O3 -std=c++0x

FILES = utils.cpp ../*.cpp ../*.h

all : ildl.mex*

MEX ?= mex
MEXFLAGS ?=

MEXFLAGS += 
MEXFLAGS += -I../ 

LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
    MEXFLAGS += -largeArrayDims
    CXXFLAGS += -DPRPACK_MEX_64
else
    CFLAGS += -DPRPACK_MEX_32
endif

MEXEXTRA ?=

MEXEXTRA += utils.cpp ../*.cpp
MEXEXTRA += CXXFLAGS="\$$CXXFLAGS $(CXXFLAGS)" LDFLAGS="\$$LDFLAGS $(CXXFLAGS)"

%.mex*: %.cpp
	$(MEX) $(MEXFLAGS) $< $(MEXEXTRA)

ildl.mex* : ildl.cpp
    
clean:
	$(RM) *.mex*
    
test: 
	matlab -nodesktop -nosplash -r test_main

debug: 
	matlab -Dgdb
.PHONY : all clean test