.PHONY: all clean test

# Project subdirectories
SRCDIR = src
INCDIR = inc
OBJDIR = obj
BINDIR = bin

# Compiler rules
# CC := icc
CC := mpicc
CFLAGS := -cc=icc -std=c11 -fPIC -I$(INCDIR) -Wall -Wextra -Werror \
	-xHOST -Ofast -ipo -no-prec-div \
	-restrict -vec-report6
# CFLAGS := -I$(INCDIR) -Wall -Wextra -Werror -Ofast -DVECTORIZE=0 \
#	-m64 -march=core-avx2 -fopt-info-vec -ftree-vectorizer-verbose=6

# Linker flags and libs
LDFLAGS := -pie
# LDLIBS := -lm
LDLIBS := -limf

# Name of final executable
TARGET := run_md

# Expanding these variables immediately with := operator
SOURCES := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(INCDIR)/*.h)
OBJECTS := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

all: $(BINDIR)/$(TARGET)

# Notes: $@ refers to the file name of the target of the rule
#        $< refers to the name of the first prerequisite
$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(LDLIBS) $(OBJECTS) -o $@

# This is a poor way to do this, actually, because an object
# file depends on all header files. I suppose we don't expect
# the headers to change very often. It would be cleaner to build
# a list of header files that each source file depends on, especially
# since this will be a relatively small project.
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(INCLUDES)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(BINDIR)/$(TARGET)
	rm -f $(OBJECTS)

test: all
	./do-sim.sh
