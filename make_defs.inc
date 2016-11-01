ifeq ($(shell uname), Darwin)
	CC := clang++
	CFLAGS := -O3 -std=c++11
	LINKFLAGS := -stdlib=libc++
	MATLAB_BIN = /Applications/Matlab/MATLAB_R2015b.app/bin/mex
else
	CC := g++
	CFLAGS := -O3 -std=gnu++0x -DNDEBUG
	LINKFLAGS := -DNDEBUG
	MATLAB_BIN = mex
endif

DEBUG := 
