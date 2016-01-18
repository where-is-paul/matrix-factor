ifeq ($(shell uname), Darwin)
	CC := clang++
	CFLAGS := -O3 -std=c++11
	LINKFLAGS := -stdlib=libc++
else
	CC := g++
	CFLAGS := -O3 -std=gnu++0x
	LINKFLAGS :=
endif

DEBUG := -w
MATLAB_BIN = /Applications/Matlab/MATLAB_R2015b.app/bin
