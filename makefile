include make_defs.inc

SRCDIR := .
BUILDDIR := build
TARGET_SYM := ldl_driver
TARGET_SKEW := skew_ldl_driver
TARBALL := matrix_factor.tar
OUTPUT := output_matrices/out*

SRC_GFLAGS = $(addprefix include/gflags/, gflags gflags_nc gflags_completions gflags_reporting)
INC_GFLAGS = $(addprefix -I, include/gflags)
INC_SYM    = $(addprefix -I, source/symmetric_matrix)
INC_SKEW   = $(addprefix -I, source/skew_symmetric_matrix)

include/gflags/%.o: include/gflags/%.cc
	$(CC) -c $(DEBUG) $(CFLAGS) $(INC_GFLAGS) $< -o $@

%.o: %.cpp
	$(CC) -c $(DEBUG) $(CFLAGS) $(INC_GFLAGS) $(INC_SYM) $(INC_SKEW) $< -o $@

$(TARGET_SYM): $(addsuffix .o, $(SRC_GFLAGS) $(TARGET_SYM))
	$(CC) $(LINKFLAGS) $? -o $@

$(TARGET_SKEW): $(addsuffix .o, $(SRC_GFLAGS) $(TARGET_SKEW))
	$(CC) $(LINKFLAGS) $? -o $@

matlab:
	cd matlab_files
	make

all: $(TARGET_SYM) $(TARGET_SKEW) matlab
	[[ -d output_matrices ]] || mkdir -p output_matrices

clean:
	$(RM) $(addsuffix .o, $(TARGET_SYM) $(TARGET_SKEW) $(SRC_GFLAGS))
	$(RM) -r $(TARGET_SYM) $(TARGET_SKEW) $(TARBALL) $(OUTPUT)

tar:
	tar cfv matrix_factor.tar ldl_driver.cpp skew_ldl_driver.cpp source

test:
	@cd matlab_files; make --no-print-directory test

