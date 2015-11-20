CC := g++
SRCDIR := $(PWD)
BUILDDIR := build
CFLAGS := -O3 -std=gnu++0x -I$(SRCDIR)/include/gflags
DEBUG := -w
TARGET := ldl_driver
TARBALL := matrix_factor.tar
OUTPUT := output_matrices/out*
 
VPATH = .:include/gflags
SOURCES_CPP := $(SRCDIR)/ldl_driver.cpp
SOURCES_CC := $(wildcard $(SRCDIR)/include/gflags/*.cc)

%.o: %.cpp
	$(CC) -c -o $@ $(DEBUG) $(CFLAGS) $<

%.o: %.cc
	$(CC) -c -o $@ $(DEBUG) $(CFLAGS) $<

$(TARGET): $(SOURCES_CPP:.cpp=.o) $(SOURCES_CC:.cc=.o)
	@mkdir -p output_matrices
	@echo " Compiling executable...";
	$(CC) -o $@ $(DEBUG) $(CFLAGS) $^
	@echo " Done.";
	@echo " Compiling mex files...";
	@cd matlab_files; make
	@echo " Done.";
	
.PHONY : clean
clean:
	@echo " Cleaning...";
	$(RM) -r $(TARGET) $(TARBALL) $(OUTPUT)
	$(RM) $(SOURCES_CPP:.cpp=.o) $(SOURCES_CC:.cc=.o)

tar:
	tar cfv matrix_factor.tar ldl_driver.cpp source

test:
	@cd matlab_files; make --no-print-directory test

-include $(DEPS)
 
.PHONY: clean
