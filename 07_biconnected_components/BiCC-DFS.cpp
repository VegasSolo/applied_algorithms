#include "stdafx.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <chrono>
#include <stack>
#include <algorithm>
#include <queue>
#ifdef _OPENMP
	#include <omp.h>
#endif
using namespace std;


/*
Vertex.h
*/

class Vertex {
public:
	int parent, rank, id, depth;
	bool discovered = false;

	Vertex();
	Vertex(const Vertex &other);
	Vertex(size_t v);

private:
	friend bool operator< (const Vertex &v1, const Vertex &v2);
	friend bool operator== (const Vertex &v1, const Vertex &v2);
	bool operator() (const Vertex &v1, const Vertex &v2) {
		return v1.id < v2.id;
	}
};

bool operator< (const Vertex &v1, const Vertex &v2) {
	return v1.id < v2.id;
}

bool operator== (const Vertex &v1, const Vertex &v2) {
	return (v1.id == v2.id);
}

/*
Vertex.cpp
*/

Vertex::Vertex()
{
	id = -1;
	parent = -1;
	rank = -1;
	depth = -1;
	discovered = false;
}

Vertex::Vertex(const Vertex & other)
{
	id = other.id;
	parent = other.id;
	rank = other.rank;
	depth = other.depth;
	discovered = other.discovered;
}

Vertex::Vertex(size_t v)
{
	id = v;
	parent = v;
	rank = 0;
	depth = -1;
	discovered = false;
}



/*
Edge.h
*/

class Edge {
public:
	int src, dest;
	string color;

	Edge();
	Edge(const Edge &other);
	Edge(size_t u, size_t v);

	Vertex &getU();
	Vertex &getV();
	pair<Vertex, Vertex> getEndpoints();

private:
	pair <Vertex, Vertex> endpoints;

	friend bool operator< (const Edge &e1, const Edge &e2);
	bool operator() (const Edge &e1, const Edge &e2) {
		if (e1.src == e2.src) {
			return e1.dest < e2.dest;
		}
		return e1.src < e2.src;
	}
};

bool operator< (const Edge &e1, const Edge &e2) {
	if (e1.src == e2.src) {
		return e1.dest < e2.dest;
	}
	return e1.src < e2.src;
}

/*
Edge.cpp
*/

Edge::Edge()
{
	src = -1;
	dest = -1;
}

Edge::Edge(const Edge & other)
{
	src = other.src;
	dest = other.dest;
	endpoints = other.endpoints;
}

Edge::Edge(size_t u, size_t v)
{
	src = u;
	dest = v;
	Vertex U(u);
	Vertex V(v);
	endpoints.first = U;
	endpoints.second = V;
}

Vertex &Edge::getU()
{
	return endpoints.first;
}

Vertex &Edge::getV()
{
	return endpoints.second;
}

pair<Vertex, Vertex> Edge::getEndpoints()
{
	return endpoints;
}


/*
Graph.h
*/

class Graph {
public:
	int numV;
	int numE;
	vector<set<Edge>> g;
	set<Vertex> vertices;

	Graph();
	Graph(int v, int e);
	Graph(const Graph &other);
	~Graph();

	void addEdge(size_t u, size_t v);

	vector<Edge> getAdjacentEdges(Vertex v);

	void Print() const;

private:
	set<Vertex>::iterator findVertex(const Vertex &vertex) const;
};

/*
Graph.cpp
*/

Graph::Graph()
{
	numV = 0;
	numE = 0;
}

Graph::Graph(int v, int e) {
	g.resize(v);
	numV = v;
	numE = e;
}

Graph::Graph(const Graph & other)
{
	g.resize(other.numV);
	numV = other.numV;
	numE = other.numE;
	vertices = other.vertices;
	for (int i = 0; i < numV; i++) {
		g[i] = other.g[i];
	}
}

Graph::~Graph() { }

void Graph::addEdge(size_t u, size_t v) {
	Edge e(u, v);
	g[u].insert(e);

	vertices.insert(e.getU());
	vertices.insert(e.getV());
}

vector<Edge> Graph::getAdjacentEdges(Vertex v)
{
	vector<Edge> edges;

	set<Vertex>::iterator it = findVertex(v);
	if (it == vertices.end())
		return edges;

	for (set<Edge>::iterator eit = g[v.id].begin(); eit != g[v.id].end(); eit++) {
		Edge e = *eit;
		edges.push_back(e);
	}

	return edges;
}

void Graph::Print() const
{
	cout << "The graph contains " << vertices.size() << " vertices:" << endl;
	set<Edge>::iterator it;
	for (int i = 0; i < vertices.size(); i++) {
		for (it = g[i].begin(); it != g[i].end(); it++) {
			Edge e = *it;
			cout << "  (" << e.src << "," << e.dest << ")" << endl;
		}
	}
}

set<Vertex>::iterator Graph::findVertex(const Vertex &vertex) const
{
	for (set<Vertex>::iterator it = vertices.begin(); it != vertices.end(); it++) {
		if (*it == vertex)
			return it;
	}
	return vertices.end();
}

/*
main.cpp
*/

Graph readFile(string fileName, vector< set<Edge> > &adj);
vector<Graph> DFS(Graph g, vector<vector<Edge>> adjE, vector<Vertex> verts);
vector<Vertex> BFS_depth(Graph g);
Graph Boruvka(Graph g, int num_threads, vector< set<Edge> > &edges);
Graph Bridges(Graph h, Graph t, vector<Vertex> verts);
void DFS_articulation(vector <int> &a, Graph g, vector<vector<Edge>> adjE);
void DFS_recurse(int u, bool visited[], int discovery[], int low[], int parent[], bool ap[], vector<vector<Edge>> adjE);
void writeResult(string fileName, vector<vector<int>> a, long long int dur, vector<int> sol);
void MakeSet(vector<Vertex> &v, int i);
void Union(vector<Vertex> &v, int x, int y);
void Link(vector<Vertex> &v, int x, int y);
int FindSet(vector<Vertex> &v, int i);

int main(int argc, char **argv)
{
	int thread_count = strtol(argv[1], NULL, 10);
	string inputFile = "graph_0.txt";
	string outputFile = "result_0.txt";

	// read graph
	vector<set<Edge>> adj;
	Graph g = readFile(inputFile, adj);
	cout << "\nTotal graph: " << endl;
	g.Print();

	// dfs to determine components
	vector <Graph> disconnectedG(g.numV);
	vector<Vertex> verts;
	set<Vertex>::iterator vit = g.vertices.begin();
	vector<vector<Edge>> adjE(g.numV);
	for (int i = 0; i < g.numV; i++) {
		if (vit != g.vertices.end())
			adjE[i] = g.getAdjacentEdges(*vit);
		verts.push_back(i);
		vit++;
	}
	disconnectedG = DFS(g, adjE, verts);


	// bfs for the purpose of setting depth correctly
	vector<vector<Vertex>> discV(disconnectedG.size());
	for (int i = 0; i < disconnectedG.size(); i++) {
		discV[i] = BFS_depth(disconnectedG[i]);
	}

	std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

	// call boruvka on each disconnected component
	vector <Graph> Tree(disconnectedG.size());
	for (int i = 0; i < disconnectedG.size(); i++) {
		Tree[i] = Boruvka(disconnectedG[i], thread_count, adj);
		cout << "\nBoruvka Tree #" << i << ":" << endl;
		Tree[i].Print();
	}

	// call bridges on each tree
	vector <Graph> component(Tree.size()); 
	for (int i = 0; i < Tree.size(); i++) {
		component[i] = Bridges(disconnectedG[i], Tree[i], verts);
	}
	
	// on each component, dfs
	vector<vector<int>> articulation(component.size());
#pragma omp parallel for num_threads(thread_count)
	for (int i = 0; i < component.size(); i++) {
		DFS_articulation(articulation[i], component[i], adjE);
	}


	std::chrono::high_resolution_clock::time_point stop_time = std::chrono::high_resolution_clock::now();
	long long int duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time - start_time).count();

	// solution creation
	vector<vector<int>> solution(1);
	for (int i = 0; i < 1; i++) {
		DFS_articulation(solution[0], g, adjE);
	}
	sort(solution[0].begin(), solution[0].end());
	solution[0].erase(unique(solution[0].begin(), solution[0].end()), solution[0].end());

	// write to file
	writeResult(outputFile, articulation, duration, solution[0]);


	return 0;
}

/* Reads an input file and returns a Graph */
Graph readFile(string fileName, vector< set<Edge> > &adj) {
	int numV = 0, numE = 0;
	int src, dest;

	ifstream infile;
	infile.open(fileName.c_str());
	if (infile.fail())
		cout << "Cannot open input file. Program will fail." << endl;

	infile >> numE;
	infile >> numV;
	Graph G(numV, numE);
	adj.resize(numE);
	Edge e;

	do {
		infile >> src;
		infile >> dest;
		e.src = src;
		e.dest = dest;

		/* IF DIRECTED, *FIXME* */
		if (src > dest)
			swap(src, dest);

		G.addEdge(src, dest);
		adj[e.src].insert(e);
		adj[e.dest].insert(e);

	} while (!infile.eof());

	infile.close();

	return G;
}

/* Runs recursive DFS on a graph to determine if it is disconnected */
vector<Graph> DFS(Graph g, vector<vector<Edge>> adjE, vector<Vertex> verts) {
	vector<Graph> components;
	Graph c(g.numV, g.numE);
	Vertex v;
	stack<Vertex> S, S2;
	bool exists = false;

	for (int i = 0; i < verts.size(); i++) {
		if (!(verts[i].discovered)) {
			S.push(verts[i]);
			while (!S.empty()) {
				v = S.top();
				S.pop();
				S2 = S;
				while (!S2.empty() && !exists) {
					if (S2.top().id == v.id)
						exists = true;
					S2.pop();
				}
				if (!v.discovered && !exists) {
					int j = 0;
					int uIndex;
					for (j = 0; j < verts.size(); j++) {
						if (v.id == verts[j].id && !verts[j].discovered) {
							verts[j].discovered = true;
							uIndex = j;
						}
					}
					for (int k = 0; k < adjE[uIndex].size(); k++) {
						for (int l = 0; l < verts.size(); l++) {
							if (adjE[v.id][k].getV() == verts[l] && !verts[l].discovered) {
								S.push(verts[l]);
								c.addEdge(verts[uIndex].id, verts[l].id);
							}
						}
					}
				}
				else if (exists)
					exists = false;
			}
			components.push_back(c);
		}
	}

	return components;
}

/* BFS for purpose of ensuring that depth is correctly set */
vector<Vertex> BFS_depth(Graph g) {
	vector<Vertex> verts;
	set<Vertex>::iterator vit = g.vertices.begin();
	vector<vector<Edge>> edges(g.numV);
	for (int i = 0; i < g.numV; i++) {
		if (vit != g.vertices.end())
			edges[i] = g.getAdjacentEdges(*vit);
		verts.push_back(i);
		if (vit != g.vertices.end())
			vit++;
	}
	Vertex current;
	
	for (int i = 0; i < verts.size(); i++) {
		verts[i].depth = 2147483647;
		verts[i].parent = -1;
	}

	queue<Vertex> Q;
	verts[0].depth = 0;
	Q.push(verts[0]);
	
	while (!Q.empty()) {
		current = Q.front();
		Q.pop();
		for (int i = 0; i < edges[current.id].size(); i++) {
			if (verts[edges[current.id][i].dest].depth == 2147483647) {
				verts[edges[current.id][i].dest].depth = current.depth + 1;
				verts[edges[current.id][i].dest].parent = current.id;
				Q.push(verts[edges[current.id][i].dest]);
			}
		}
	}
	return verts;
}

/* Uses Boruvka algorithm to find a Tree from the Graph and returns it */
Graph Boruvka(Graph g, int thread_count, vector< set<Edge> > &edges) {
	Graph Tree(g.numV, g.numE);

	Edge e;
	vector <Vertex> v(g.numV);
	vector <Edge> activeEdges(g.numV);
	int numAE = g.numV;
	for (int i = 0; i < g.numV; i++) {
		MakeSet(v, i);
		if (edges[i].begin() != edges[i].end()) {
			activeEdges[i] = *edges[i].begin();
			edges[i].erase(activeEdges[i]);
		}
		else
			numAE--;
		if (edges[i].begin() == edges[i].end())
			numAE--;
	}
	int numComponents = g.numV;
	

	set <Edge> mstE;
	set <Edge> mstEdges;

	while (numComponents > 1 && numAE > 0) {
#pragma omp parallel for num_threads(thread_count)
		for (int i = 0; i < g.numV; i++) {
			Edge foundE = activeEdges[i];
			if (FindSet(v, foundE.src) != FindSet(v, foundE.dest)) {
				mstE.insert(foundE);
			}
			else {
				bool found = false;
				while (!edges[i].empty() && !found) {
					edges[i].erase(activeEdges[i]);
					if (!edges[i].empty()) {
						activeEdges[i] = *edges[i].begin();
					}
					else
						numAE--;
					foundE = activeEdges[i];
					if (FindSet(v, foundE.src) != FindSet(v, foundE.dest)) {
						mstE.insert(foundE);
						found = true;
					}
				}
			}
		}
		for (set<Edge>::iterator it = mstE.begin(); it != mstE.end(); it++) {
			e = *it;
			if (FindSet(v, e.src) != FindSet(v, e.dest)) {
				Union(v, e.src, e.dest);
				mstEdges.insert(e);
				numComponents--;
				if (!edges[e.src].empty()) {
					activeEdges[e.src] = *edges[e.src].begin();
					edges[e.src].erase(e);
				}
				else
					numAE--;
				if (!edges.at(e.dest).empty()) {
					activeEdges.at(e.dest) = *edges.at(e.dest).begin();
					edges.at(e.dest).erase(e);
				}
				else
					numAE--;
			}
			else {
				mstE.erase(e);
			}
		}
	}
	
	for (std::set<Edge>::iterator it = mstEdges.begin(); it != mstEdges.end(); it++) {
		Edge e = *it;
		Tree.addEdge(e.src, e.dest);
	}
	return Tree;
}

/* Returns a graph that was a tree with bridges removed */
Graph Bridges(Graph h, Graph t, vector<Vertex> verts) {
	Graph B(h.numV, h.numE);
	Graph GminusB(h.numV, h.numE);
	Graph GminusT(h.numV, h.numE);
	bool isInT = false;
	bool isInB = false;

	for (int i = 0; i < h.numV; i++)
	{
		for (std::set<Edge>::iterator it = h.g[i].begin(); it != h.g[i].end(); it++)
		{
			Edge e = *it;
			isInT = false;
			for (int j = 0; j < t.numV; j++)
			{
				for (std::set<Edge>::iterator ite = t.g[j].begin(); ite != t.g[j].end(); ite++)
				{
					Edge d = *ite;
					if (e.getU() == d.getU() && e.getV() == d.getV())
					{
						isInT = true;
						break;
					}
				}
				if (isInT == true)
				{
					break;
				}
			}
			if (isInT == false)
			{
				GminusT.addEdge(e.src, e.dest);
			}
		}
	}
	for (int i = 0; i < GminusT.numV; i++)
	{
		for (std::set<Edge>::iterator it = GminusT.g[i].begin(); it != GminusT.g[i].end(); it++)
		{
			Edge e = *it;
			Vertex u = e.getU();
			Vertex v = e.getV();
			Vertex U = verts[u.id];
			Vertex V = verts[v.id];
			while (U.depth != V.depth)
			{
				h.g[u.id].erase(e);
				h.vertices.erase(e.getU());
				h.vertices.erase(e.getV());
				e.color = "green";
				if (U.depth < V.depth)
				{
					U.id = U.parent;
				}
				else
				{
					V.id = V.parent;
				}
				h.g[u.id].insert(e);
				h.vertices.insert(U);
				h.vertices.insert(V);
			}
			while (U.parent != V.parent && U.id != V.id)
			{
				h.g[u.id].erase(e);
				h.vertices.erase(e.getU());
				h.vertices.erase(e.getV());
				e.color = "green";
				U.id = U.parent;
				V.id = V.parent;
				h.g[u.id].insert(e);
				h.vertices.insert(U);
				h.vertices.insert(V);
			}
		}
	}
	for (int i = 0; i < h.numV; i++)
	{
		for (std::set<Edge>::iterator it = h.g[i].begin(); it != h.g[i].end(); it++)
		{
			Edge e = *it;
			if (e.color == "black")
			{
				B.addEdge(e.src, e.dest);
			}
		}
	}
	for (int i = 0; i < h.numV; i++)
	{
		for (std::set<Edge>::iterator it = h.g[i].begin(); it != h.g[i].end(); it++)
		{
			Edge e = *it;
			isInB = false;
			for (int j = 0; j < B.numV; j++)
			{
				for (std::set<Edge>::iterator ite = B.g[j].begin(); ite != B.g[j].end(); ite++)
				{
					Edge d = *ite;
					if (e.getU() == d.getU() && e.getV() == d.getV())
					{
						isInB = true;
						break;
					}
				}
				if (isInB == true)
				{
					break;
				}
			}
			if (isInB == false)
			{
				GminusB.addEdge(e.src, e.dest);
			}
		}
	}
	return GminusB;
}

/* Runs recursive DFS to return the articulation points */
void DFS_articulation(vector <int> &a, Graph g, vector<vector<Edge>> adjE) {
	bool *visited = new bool[g.numV];
	int *discovery = new int[g.numV];
	int *low = new int[g.numV];
	int *parent = new int[g.numV];
	bool *ap = new bool[g.numV];

	for (int i = 0; i < g.numV; i++) {
		parent[i] = -1;
		visited[i] = false;
		ap[i] = false;
	}
	for (int i = 0; i < g.numV; i++) {
		if (visited[i] == false)
			DFS_recurse(i, visited, discovery, low, parent, ap, adjE);
	}
	for (int i = 0; i < g.numV; i++) {
		if (ap[i] == true)
			a.push_back(i);
	}
}

/* Recursive utility for DFS to find articulation points */
void DFS_recurse(int u, bool visited[], int discovery[], int low[], int parent[], bool ap[], vector<vector<Edge>> adjE) {
	static int time = 0;
	int children = 0;
	visited[u] = true;
	discovery[u] = low[u] = ++time;
	for (int i = 0; i < adjE[u].size(); i++) {
		int v = adjE[u][i].dest;
		if (!visited[v]) {
			children++;
			parent[v] = u;
			DFS_recurse(v, visited, discovery, low, parent, ap, adjE);
			low[u] = min(low[u], low[v]);
			if (parent[u] == -1 && children > 1)
				ap[u] = true;
			if (parent[u] != -1 && low[v] >= discovery[u])
				ap[u] = true;
			else if (v != parent[u])
				low[u] = min(low[u], discovery[v]);
		}
	}
}

/* Writes the results of timing and articulation points to file */
void writeResult(string fileName, vector<vector<int>> a, long long int dur, vector<int> sol) {
	ofstream outfile;
	outfile.open(fileName.c_str());
	if (outfile.fail())
		cout << "Cannot open output file. Program will fail." << endl;

	outfile << "The program completed in " << (dur / 1E9) << " seconds" << endl;
	outfile << "It found these articulation points: " << endl;
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < a[i].size(); j++) {
			outfile << a[i][j] << " ";
		}
		outfile << endl;
	}
	outfile << "Expected Articulation Points: " << endl;
	for (int i = 0; i < sol.size(); i++) {
		outfile << sol[i] << " ";
	}
	outfile << "\nEnd" << endl;
}

/*
Disjoint Data Set Functions given in book
*/
void MakeSet(vector<Vertex> &v, int i) {
	v.at(i).id = i;
	v.at(i).parent = i;
	v.at(i).rank = 0;
}
void Union(vector<Vertex> &v, int x, int y) {
	Link(v, FindSet(v, x), FindSet(v, y));
}
void Link(vector<Vertex> &v, int x, int y) {
	if (v.at(x).rank > v.at(y).rank)
		v.at(y).parent = x;
	else
		v.at(x).parent = y;
	if (v.at(x).rank == v.at(y).rank)
		v.at(y).rank = v.at(y).rank + 1;
}
int FindSet(vector<Vertex> &v, int i) {
	if (v.at(i).parent != v.at(i).id) {
		v.at(i).parent = FindSet(v, v.at(i).parent);
	}
	return v.at(i).parent;
}