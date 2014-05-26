CC := g++
SRCDIR := .
BUILDDIR := build
CFLAGS := -O3 -flto -std=gnu++0x
DEBUG := -w
TARGET_SYM := ldl_driver
TARGET_SKEW := skew_ldl_driver
TARBALL := matrix_factor.tar
OUTPUT := output_matrices/out*
 
SOURCES_SYM := ./ldl_driver.cpp ./include/gflags/*
SOURCES_SKEW := ./skew_ldl_driver.cpp ./include/gflags/*

all: 
	@mkdir -p output_matrices
	@echo " Compiling symmetric executable..."; $(CC) $^ $(DEBUG) $(CFLAGS) $(SOURCES_SYM) -o $(TARGET_SYM)
	@echo " Done.";
	@echo " Compiling skew-symmetric executable..."; $(CC) $^ $(DEBUG) $(CFLAGS) $(SOURCES_SKEW) -o $(TARGET_SKEW)
	@echo " Done.";
	@echo " Compiling mex files...";
	@cd matlab_files; make > /dev/null
	@echo " Done.";
	
.PHONY : clean
clean:
	@echo " Cleaning..."; $(RM) -r $(TARGET_SYM) $(TARGET_SKEW) $(TARBALL) $(OUTPUT); 

tar:
	tar cfv matrix_factor.tar ldl_driver.cpp source

test:
	@cd matlab_files; make --no-print-directory test

-include $(DEPS)
 
.PHONY: clean
