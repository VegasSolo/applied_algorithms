/*
	Vertex.h
*/

class Vertex {
public:
	int parent, rank, id, depth;
private:
	friend bool operator< (const Vertex &v1, const Vertex &v2);
	bool operator() (const Vertex &v1, const Vertex &v2) {
		return v1.id < v2.id;
	}
};

bool operator< (const Vertex &v1, const Vertex &v2) {
	return v1.id < v2.id;
}
