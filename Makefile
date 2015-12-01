# See LICENSE.txt for license details.

PAR_FLAG = -fopenmp

ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG = -openmp
endif

CXX_FLAGS += -std=c++11 -O3 -Wall $(PAR_FLAG)

ifneq (,$(findstring sunCC,$(CXX)))
	CXX_FLAGS = -std=c++11 -xO3 -m64 -xtarget=native -xopenmp
endif

KERNELS = bc bfs cc pr sssp tc
SUITE = $(KERNELS) converter

.PHONY: all
all: $(SUITE)

% : src/%.cc src/*.h
	$(CXX) $(CXX_FLAGS) $< -o $@$

# Testing
include test/test.mk

.PHONY: clean
clean:
	rm -f $(SUITE) test/out/*
