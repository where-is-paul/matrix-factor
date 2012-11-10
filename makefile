# Generic Makefile for compiling a simple executable.
 
CC := g++
SRCDIR := .
BUILDDIR := build
CFLAGS := -O3 -funroll-loops -std=c++0x
DEBUG := -Wextra -Wall
TARGET := ldl_driver
TARBALL := matrix_factor.tar
OUTPUT := output_matrices/out*
 
SRCEXT := cpp
SOURCES := ./ldl_driver.cpp
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
DEPS := $(OBJECTS:.o=.deps)
	
$(TARGET): $(OBJECTS)
	@echo " Linking..."; $(CC) $^ $(DEBUG) -o $(TARGET)
	@cd matlab_files; make --no-print-directory
 
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p output_matrices
	@mkdir -p $(BUILDDIR)
	@echo " CC $<"; $(CC) $(DEBUG) $(CFLAGS) -MD -MF $(@:.o=.deps) -c -o $@ $<
	
clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET) $(TARBALL) $(OUTPUT); 

tar:
	tar cfv matrix_factor.tar ldl_driver.cpp source

test:
	@cd matlab_files; make --no-print-directory test

-include $(DEPS)
 
.PHONY: clean
