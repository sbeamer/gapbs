# See LICENSE.txt for license details.

CXX_FLAGS += -std=c++11 -O3 -fopenmp

SUITE = bc bfs cc pr sssp tc converter

.PHONY: all
all: $(SUITE)

% : %.cc *.h
	$(CXX) $(CXX_FLAGS) $< -o $@$

.PHONY: clean
clean:
	rm -f $(SUITE)
