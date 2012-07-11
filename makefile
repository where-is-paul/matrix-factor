# Generic Makefile for compiling a simple executable.
 
CC := g++
SRCDIR := .
BUILDDIR := build
CFLAGS := -Wall -O3 -funroll-loops -std=c++0x
TARGET := ldl_driver
 
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -name *.$(SRCEXT))
ALL := $(shell find $(SRCDIR) -name "*.h" -or -name "*.cpp")
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
DEPS := $(OBJECTS:.o=.deps)
 
$(TARGET): $(OBJECTS)
	@echo " Linking..."; $(CC) $^ -o $(TARGET)
 
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " CC $<"; $(CC) $(CFLAGS) -MD -MF $(@:.o=.deps) -c -o $@ $<
 
clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET)

tar:
	tar cfv matrix_factor.tar $(ALL)

-include $(DEPS)
 
.PHONY: clean
