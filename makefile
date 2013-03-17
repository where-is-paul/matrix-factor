CC := g++
SRCDIR := .
BUILDDIR := build
CFLAGS := -O3 -std=c++0x
DEBUG := -Wextra -Wall -pg
TARGET := ldl_driver
TARBALL := matrix_factor.tar
OUTPUT := output_matrices/out*
 
SOURCES := ./ldl_driver.cpp

$(TARGET): 
	@mkdir -p output_matrices
	@echo " Compiling executable..."; $(CC) $^ $(DEBUG) $(CFLAGS) $(SOURCES) -o $(TARGET)
	@echo " Done.";
	@echo " Compiling mex files...";
	@cd matlab_files; make > /dev/null
	@echo " Done.";
	
.PHONY : clean
clean:
	@echo " Cleaning..."; $(RM) -r $(TARGET) $(TARBALL) $(OUTPUT); 

tar:
	tar cfv matrix_factor.tar ldl_driver.cpp source

test:
	@cd matlab_files; make --no-print-directory test

-include $(DEPS)
 
.PHONY: clean
