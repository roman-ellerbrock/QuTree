#include "Vertex2D.h"

double distance(const Vertex& v1, const Vertex& v2) {
	return sqrt(pow(v1.x - v2.x, 2) + pow(v1.y - v2.y, 2));
}


