/*
	Edge.h
*/

class Edge {
public:
	int src, dest;
    string color;
	
private:
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