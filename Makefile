CXX ?= g++
MPICXX ?= mpic++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -fno-exceptions -DOMPI_SKIP_MPICXX
GRID_N ?= 4000

all: seq parallel

seq: Lab2_sequence.cpp
	$(CXX) $(CXXFLAGS) -DN=$(GRID_N) Lab2_sequence.cpp -o Lab2_sequence

parallel: Lab2_parallel.cpp
	$(MPICXX) $(CXXFLAGS) -DN=$(GRID_N) Lab2_parallel.cpp -o Lab2_parallel

clean:
	rm -f Lab2_sequence Lab2_parallel

.PHONY: all seq parallel clean
