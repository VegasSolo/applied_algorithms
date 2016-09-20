/*
Author:		Dillon VanBuskirk
Date:		6-12-2016
Program:	Boruvka's Algorithm [Parallel]
Purpose:	Read graph from file and report the sum of the edge weights for the
			Minimum Spanning Tree and all of the edges in lexicographical order.
			This will be done in parallel in order to prove it is in NC.
*/

//#include "stdafx.h" // VisualStudio
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#ifdef _OPENMP
	#include <omp.h>
#endif
using namespace std;


struct Edge {
	int weight, src, dest;

	friend bool operator< (const Edge &e1, const Edge &e2);
	bool operator() (const Edge &e1, const Edge &e2) {
		if (e1.weight == e2.weight) {
			if (e1.src == e2.src) {
				return e1.dest < e2.dest;
			}
			return e1.src < e2.src;
		}
		return e1.weight < e2.weight;
	}
};

bool operator< (const Edge &e1, const Edge &e2) {
	if (e1.weight == e2.weight) {
		if (e1.src == e2.src) {
			return e1.dest < e2.dest;
		}
		return e1.src < e2.src;
	}
	return e1.weight < e2.weight;
}

struct Vertex {
	int parent, rank, id; // do I want a pointer to the edge set? Nope.
};


void boruvka(int thread_count);
void MakeSet(vector<Vertex> &v, int i);
void Union(vector<Vertex> &v, int x, int y);
void Link(vector<Vertex> &v, int x, int y);
int FindSet(vector<Vertex> &v, int i);
void writeMST(int mst, set<Edge> mstEdges);


int main (int argc, char* argv[]) {
	int thread_count = strtol(argv[1], NULL, 10);

	boruvka(thread_count);

    return 0;
}

void boruvka(int thread_count) {
	// declarations
	string const f = "graph1.txt";
	int numV = 0, numE = 0;

	// read file
	ifstream infile;
	infile.open(f.c_str());
	if (infile.fail())
		cout << "Cannot open input file. Program will fail." << endl;

	int index = 0;
	Edge e;
	infile >> numV;
	infile >> numE;

	vector < set<Edge> > edges(numE); // This is the master container for all the edges. The vector at index i contains the set of Edges for vertex i
	//vector < set<Edge> > edgesBackup(numE); // Backup set of edges which will not be dynamically deleted as edges are used. This was never needed in my code.

	do {
		infile >> e.src;
		infile >> e.dest;
		infile >> e.weight;

		edges.at(e.src).insert(e); // This inserts the edge <src, dest, weight> into the set at the vector position that corresponds to vertex of id src
		edges.at(e.dest).insert(e); // This inserts the edge <src, dest, weight> into the set at the vector position that corresponds to vertex of id dest
		//edgesBackup.at(e.src).insert(e); // These two are the backup in order to maintain 
		//edgesBackup.at(e.dest).insert(e); // a complete list of edges for each vertex

		index++;

	} while (!infile.eof());

	//cout << "Finished reading. Closed the file & calling MakeSet & populating activeEdges" << endl;
	infile.close();

	// MakeSet on all vertices and populate the lightweight edge set, activeEdges
	vector <Vertex> v(numV); // This is a vector of vertices that contains the rank, parent, and id of each vertex at index i of the vector
	vector<Edge> activeEdges(numV); // This is a vector of edges that are identified as the lightweight edge out of each index i of the vector (vertex i)
	int numAE = numV; // This keeps track of the count of activeEdges. It is decremented when there are no more edges from master set that can be added to activeEdges
	for (int i = 0; i < numV; i++) {
		MakeSet(v, i);
		activeEdges.at(i) = *edges.at(i).begin(); // activeEdges of vertex i = the beginning of the set that is in the ith index of the vector
		edges.at(i).erase(activeEdges.at(i)); // erases the element of the set that we just placed in activeEdges (the chepaest edge in the set) of the ith index of the vector
	}
	int numComponents = numV; // This is the number of components. It is decremented after each Union

	// create the MST container
	set<Edge> mstE; // This is the set of edges that is identified as lightweight edges out of components. It is dynamic and is emptied at the end of each iteration of the while loop
	set<Edge> mstEdges; // This is the set of edges that have been added to the MST. This set is not changed once edges have been added to it. The additions occur after each Union
	int mst = 0; // This is the running sum of the weights of each edge that is added to mstEdges. 

	// while there is more than one component and while there exists at least one lightweight edge out of a component
	while (numComponents > 1 && numAE > 0) {
		// parallel OMP for
		// for each component
#pragma omp parallel for
		for (int i = 0; i < numV; i++) {
#ifdef _OPENMP 
			int my_rank = omp_get_thread_num();
			int thread_count = omp_get_num_threads();
#else
			int my_rank = 0;
			int thread_count = 1;
#endif
			cout << "Hi, I'm thread " << my_rank << " of " << thread_count << endl;

			// find lightweight edge <u,v>
			Edge lightE = activeEdges.at(i);
			cout << "Found lightweight edge for " << i << "th vertex as (" << lightE.src << "," << lightE.dest << "," << lightE.weight << ")" << endl;
			// if FindSet(e.src) == FindSet(e.dest) do nothing. this is the same component. else add the lightweight edge to MST
			if (FindSet(v, lightE.src) != FindSet(v, lightE.dest)) { 
				// OMP lock --- not needed because of the way I set this up. Yay!
				mstE.insert(lightE);
				// OMP unlock
				cout << "That edge's vertices are not in the same component. It was added to the MST." << endl;
			}
			else {
				cout << "That edge's vertices are in the same component. Removing that edge." << endl;
				// OMP lock ? --- not needed because of the way I set this up. Yay!
				bool found = false;
				while (!edges.at(i).empty() && !found) {
					edges.at(i).erase(activeEdges.at(i));
					if (!edges.at(i).empty()) {
						activeEdges.at(i) = *edges.at(i).begin();
						cout << "Update " << i << " activeEdges to (" << activeEdges.at(i).src << "," << activeEdges.at(i).dest << "," << activeEdges.at(i).weight << ")" << endl;
					}
					else
						numAE--;
					lightE = activeEdges.at(i);
					if (FindSet(v, lightE.src) != FindSet(v, lightE.dest)) {
						mstE.insert(lightE);
						cout << "That edge's vertices are not in the same component. It was added to the MST." << endl;
						found = true;
					}
				}
				// OMP unlock
			}
		}
		// for each edge in mstE (the set of lightweight edges out of every vertex)
		for (std::set<Edge>::iterator it = mstE.begin(); it != mstE.end(); it++) {
			e = *it;
			// if the src and destination are not in the same component. If they are not, then Union(src,dest), add to mstEdges, decrement numC, add weight to mst
			if (FindSet(v, e.src) != FindSet(v, e.dest)) {
				cout << "From the MST, Union " << e.src << " and " << e.dest << " of " << e.weight << endl;
				Union(v, e.src, e.dest);
				mst += e.weight;
				mstEdges.insert(e);
				numComponents--;
				// Update edges here. 
				// If the set at vertex i is not empty, then add the next edge (which is the lightweight edge for that vertex) to activeEdges.
				if (!edges.at(e.src).empty()) {
					activeEdges.at(e.src) = *edges.at(e.src).begin();
					edges.at(e.src).erase(e);
					cout << "Update " << e.src << " activeEdges to (" << activeEdges.at(e.src).src << "," << activeEdges.at(e.src).dest << "," << activeEdges.at(e.src).weight << ")" << endl;
				}
				// Else, decrement AE because that vertex doesn't have anymore edges
				else
					numAE--;
				if (!edges.at(e.dest).empty()) {
					activeEdges.at(e.dest) = *edges.at(e.dest).begin();
					edges.at(e.dest).erase(e);
					cout << "Update " << e.dest << " activeEdges to (" << activeEdges.at(e.dest).src << "," << activeEdges.at(e.dest).dest << "," << activeEdges.at(e.dest).weight << ")" << endl;
				}
				else
					numAE--;
			}
			else {
				// if the edge that is in mstE has two vertices that are in the same component, then remove that edge from mstE and continue
				mstE.erase(e);
				it = mstE.begin();
			}
		}
	}

	cout << "\nThe final weight of the MST is " << mst << endl;
	// write to file
	writeMST(mst, mstEdges);
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

/*
Write sum of MST and list of MST edges in lexicographical order
*/
void writeMST(int mst, set<Edge> mstEdges) {
	string const outputFile = "mst.txt";

	ofstream outfile;
	outfile.open(outputFile.c_str());
	if (outfile.fail())
		cout << "Cannot open output file. Program will fail." << endl;

	outfile << "The sum of MST is " << mst << endl;
	outfile << "The list of edges in the MST is: " << endl;

	for (std::set<Edge>::iterator it = mstEdges.begin(); it != mstEdges.end(); it++) {
		Edge e = *it;
		outfile << "(" << e.src << "," << e.dest << "," << e.weight << ")" << endl;
	}
	outfile.close();
}
