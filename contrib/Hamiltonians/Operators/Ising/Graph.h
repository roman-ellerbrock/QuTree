#pragma once
#include "Vertex2D.h"
#include <vector>
#include <iostream>
#include <assert.h>

using namespace std;

class Edge {
	public:
	Edge(size_t src_, size_t dest_, double weight_)
		:src(src_), dest(dest_), weight(weight_) {
	}
	~Edge() = default;
	size_t src, dest;
	double weight;

	void print()const {
		cout << "( " << src << " -> " << dest << " | " << weight << " )\n";
	}
};

class Graph {
public:
	Graph() = default;
	~Graph() = default;

	void push_back(const Vertex& v) { vertices.push_back(v); }
	void push_back(const Edge& e) { edges.push_back(e); }

	Vertex vertex(const size_t i)const {
		assert(i < vertices.size());
		return vertices[i];
	}

	Edge edge(const size_t i)const {
		assert(i < edges.size());
		return edges[i];
	}

	void all_distance_edges() {
		for (size_t i = 0; i < vertices.size(); ++i) {
			for (size_t j = 0; j < vertices.size(); ++j) {
				edges.push_back(Edge(i, j, distance(vertex(i), vertex(j))));
			}
		}
	}

	void print()const {
//		for (const Vertex& v : vertices) { v.print(); }
	
		cout << "#!/bin/gnuplot" << endl;
		for (size_t nu = 0; nu < vertices.size(); ++nu) {
			const Vertex v = vertices[nu];
			cout << "set label \"" << nu << "\" at " << v.x << ", " << v.y << endl;
		}

		cout << "set xr [0:1]" << endl;
		cout << "set yr [0:1]" << endl;
		cout << "f(x)=0" << endl;
		cout << "plot f(x)" << endl;
		cout << "pause -1" << endl << endl;

		for (const Edge& e : edges) { e.print(); }
	}

	size_t nEdges()const { return edges.size(); }
	size_t nVertices()const { return vertices.size(); }

private:
	vector<Vertex> vertices;
	vector<Edge> edges;
};

