#pragma once
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <random>

using namespace std;

class Vertex {
public:
	Vertex(mt19937& gen) {
		uniform_real_distribution<double> dist(0., 1.);
		x = dist(gen);
		y = dist(gen);
	}
	Vertex(double x_, double y_): x(x_), y(y_) {}
	double x, y;

	void print()const {
		cout << "( " << x << ", " << y << " )\n";
	}
};


double distance(const Vertex& v1, const Vertex& v2);

