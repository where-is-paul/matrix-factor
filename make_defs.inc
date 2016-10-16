ifeq ($(shell uname), Darwin)
	CC := clang++
	CFLAGS := -O3 -std=c++11
	LINKFLAGS := -stdlib=libc++ -lgfortran
	MATLAB_BIN = /Applications/Matlab/MATLAB_R2015b.app/bin/mex
	MEX_EXT = $(shell $(MATLAB_BIN)/mexext)
else
	CC := g++
	CFLAGS := -O3 -std=gnu++0x -DNDEBUG
	LINKFLAGS := -lgfortran -O3 -DNDEBUG
	MATLAB_BIN = mex
	MEX_EXT = 
endif

DEBUG := -w
