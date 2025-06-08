# -----------------------------------------------------------------------
# Author: Bojun Zhang
# Date  : 2025-04-15 19:40:07
# File  : Makefile
# -----------------------------------------------------------------------
# Compiler and flags
FC      = gfortran
FFLAGS  = -O2 

# Executable name
EXEC = advection_solver

# Source and object files
SRC  = $(wildcard *.f90)
OBJ  = $(patsubst %.f90, %.o, $(SRC))

# Default target
all: $(EXEC)

# Linking
$(EXEC): $(OBJ)
	$(FC) -o $@ $^

# Compilation rule
%.o: %.f90
	$(FC) $(FFLAGS) -c $^

# Clean rule
clean:
	rm -f $(EXEC) $(OBJ) fort.*  *.txt *.dat

# Run the program (optional helper)
run: $(EXEC)
	./$(EXEC)

echo:
	@echo $(SRC)
	@echo $(OBJ)
