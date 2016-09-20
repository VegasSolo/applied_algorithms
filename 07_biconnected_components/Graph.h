/*
	Graph.h
*/

#include "Vertex.h"
#include "Edge.h"

#include <vector>
#include <set>

class Graph {
public:
	Graph();
	Graph(const Graph &other);
	~Graph();
	
	void addEdge(size_t u, size_t v);
	
	
	
private:
	int numV;
	int numE;

	vector < set <Edge> > g;
	set <Vertex> vertices;
}