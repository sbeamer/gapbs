# See LICENSE.txt for license details.

PAR_FLAG = -fopenmp

ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG = -openmp
endif

CXX_FLAGS += -std=c++11 -O3 -Wall $(PAR_FLAG)

KERNELS = bc bfs cc pr sssp tc
SUITE = $(KERNELS) converter

.PHONY: all
all: $(SUITE)

% : %.cc *.h
	$(CXX) $(CXX_FLAGS) $< -o $@$

# Testing
include test/test.mk

.PHONY: clean
clean:
	rm -f $(SUITE) test/out/*
