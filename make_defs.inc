ifeq ($(shell uname), Darwin)
	CC := clang++
	CFLAGS := -O3 -std=c++11
	LINKFLAGS := -stdlib=libc++ 
	MATLAB_BIN = /Applications/Matlab/MATLAB_R2015b.app/bin/mex
	MEX_EXT = $(shell $(MATLAB_BIN)/mexext)
else
	CC := g++
	CFLAGS := -std=gnu++0x
	LINKFLAGS := 
	MATLAB_BIN = mex.bat
	MEX_EXT = 
endif

DEBUG := -w
