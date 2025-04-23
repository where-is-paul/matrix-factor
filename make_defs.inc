ifeq ($(shell uname), Darwin)
	CC := clang++
	CFLAGS := -O3 -std=c++11 -fopenmp
	LINKFLAGS := -stdlib=libc++ -fopenmp
	MATLAB_BIN = /Applications/Matlab/MATLAB_R2015b.app/bin/mex
	MEX_EXT = $(shell $(MATLAB_BIN)/mexext)
else
	CC := g++
	CFLAGS := -O3 -std=gnu++0x -DSYM_ILDL_DEBUG -fopenmp
	LINKFLAGS := -O3 -DSYM_ILDL_DEBUG -fopenmp
	MATLAB_BIN = mex
	MEX_EXT = 
endif

DEBUG := -w
