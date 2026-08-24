CXX = g++
CXXFLAGS = -std=c++23 -O3 -march=native -DNDEBUG -fopenmp -flto -I../seqan3/include -I../xxHash
LDFLAGS = -lz -lbz2
SRC = ../code

all: function_tests generate_clusters generate_layout

function_tests: ${SRC}/fast_construct_test.cpp ${SRC}/fast_construct_hashing.cpp ${SRC}/fast_construct_graph.cpp ${SRC}/fast_construct_bin.cpp
	${CXX} ${CXXFLAGS} $^ -o $@ $(LDFLAGS)

generate_clusters: ${SRC}/fast_construct_graph.cpp ${SRC}/test_clustering.cpp ${SRC}/fast_construct_hashing.cpp ${SRC}/generate_clusters.cpp ${SRC}/fast_construct_bin.cpp
	${CXX} ${CXXFLAGS} $^ -o $@ ${LDFLAGS}

generate_layout: ${SRC}/fast_construct_graph.cpp ${SRC}/fast_construct_hashing.cpp ${SRC}/generate_layout.cpp ${SRC}/fast_construct_bin.cpp
	${CXX} ${CXXFLAGS} $^ -o $@ ${LDFLAGS}

clean:
	rm -f function_tests generate_clusters generate_layout

.PHONY: all clean
