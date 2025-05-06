# Define variables
CC = mpic++                       # MPI compiler
CFLAGS = -g -pg -fopenmp          # Compile with debugging symbols and profiling enabled
EXEC = mpi.exe                    # Output executable file
SRC = V1.cpp                      # Source file
PROFILE_FILE = gmon.out           # Profiling data file
TXT_FILE = analysis.txt           # Text file for gprof output
GPROF2DOT = python3 gprof2dot.py

all:
	@echo "Compiling the program..."
	@$(CC) $(CFLAGS) $(SRC) -o $(EXEC)
	@echo "Running the program..."
	@/usr/bin/time -f "\
Elapsed Time (Real): %E\n\
User CPU Time: %U seconds\n\
System CPU Time: %S seconds\n\
CPU Usage: %P\n\
Max Resident Memory: %M KB\n\
Avg Text Size: %X KB\n\
Avg Data Size: %D KB\n\
Inputs: %I\n\
Outputs: %O\n\
Major Page Faults: %F\n\
Minor Page Faults: %R\n\
Swaps: %W\n" \
		mpirun --hostfile hosts.txt -np 1 ./$(EXEC) -fopenmp

# Profile target: compile, run with profiling, and generate call graph
profile:
	@echo "Compiling the program with profiling..."
	@$(CC) $(CFLAGS) $(SRC) -o $(EXEC)
	@echo "Running the program with profiling..."
	@/usr/bin/time -f "\
User CPU Time: %U seconds\n\
System CPU Time: %S seconds\n\
Elapsed Time (Real): %E\n\
CPU Usage: %P\n\
Max Resident Memory: %M KB\n\
Avg Text Size: %X KB\n\
Avg Data Size: %D KB\n\
Inputs: %I\n\
Outputs: %O\n\
Major Page Faults: %F\n\
Minor Page Faults: %R\n\
Swaps: %W\n" \
	mpirun --hostfile hosts.txt -np 1 ./$(EXEC) -fopenmp
	@echo "Generating profiling data..."
	@echo "Profiling data generated in $(PROFILE_FILE)"
	@# Generate gprof analysis and callgraph
	@gprof $(EXEC) $(PROFILE_FILE) > $(TXT_FILE)
	@$(GPROF2DOT) -s -n 2 -e 2 < $(TXT_FILE) > callgraph.dot
	@dot -Tpng callgraph.dot -o callgraph.png

# Clean up generated files
clean:
	@echo "Cleaning up..."
	@rm -f $(EXEC) $(PROFILE_FILE) $(TXT_FILE) callgraph.dot callgraph.png gmon.out

# Declare phony targets
.PHONY: all profile clean
