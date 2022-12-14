#include "Graph.h"
#include <iostream>
#include <random>
#include <chrono>


int main() {
	size_t N = 5;
	using namespace std::chrono;
	auto seed = system_clock::now().time_since_epoch().count();
	mt19937 gen(seed);

	// Build a (random) Graph
	Graph graph;
	for (size_t i = 0; i < N; ++i) {
		Vertex v(gen);
		graph.push_back(v);
	}
	graph.all_distance_edges();

	graph.print();

	return 0;
}


