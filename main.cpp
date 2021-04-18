#include <iostream>
#include <time.h>
#include <vector>
#include <fstream>
#include <queue>

using std::vector;
using std::priority_queue;
using std::cout;
using std::cin;
using std::endl;

struct MST {
public:
    float obj_val; 
    vector<int> edges;

    MST() {
        obj_val = 0;
    }

    void add_edge(int _edge_index) {
        edges.push_back(_edge_index);
    }
    
    void set_objective(float _new_obj_val) {
        obj_val = _new_obj_val;
    }

    void clear() {
        obj_val = 0;
        edges.clear();
    }
};

class Edge {
public:
    int edge_index;
    int i;
    int j;
    float H;

    Edge(int _edge_index, int _i, int _j, float _H){
        edge_index = _edge_index;
        i = _i;
        j = _j;
        H = _H;
    }

    void update_Edge(float _new_H) {
        H = _new_H;
    }
};

class Vertex {
public:
    int i;
    float d_i;

    Vertex(int _i, float _d_i) {
        i = _i;
        d_i = _d_i;
    }

    void update_Vertex(float _new_d_i) {
        d_i = _new_d_i;
    }

    float getD() const { return d_i; }
};

class Comparator
{
public:
    int operator() (const Vertex& v, const Vertex& u)
    {
        return v.getD() > u.getD();
    }
};


void mst_prims(int n, int m, float C, vector<int>& neighbors, vector<float>& h, vector<int>& indexes, MST& t_star) {
    /*
     * Subroutine for finding MST(k)
     * 
     * Receive current graph topology and costs 
     * Receive MST address to update
     * Initialize --> heap, distance labels, pred, T^{*} (empty)
     * While T^{*}.edges is not has size less than n-1, pop next node, update heap for adjacent 
     * Returns T^{*} as a tree struct 
     */
    
    // pre initialization
    // t_star is the MST (indexes of edges as per lists neighbors, h, c, tau)
    // t_star has attributes edges and obj_val
    // other values - , 
    t_star->clear(); 
    
        
    // initialize heap of vertices 
    // node 0 has d(0) = 0
    // everyone else has d(j) = C,
    priority_queue <Vertex, vector<Vertex>, Comparator> heap;

    for (int i = 0; i < n; i++) {
        if (i == 0) {
            heap.push(Vertex(i, 0));
        }
        else {
            heap.push(Vertex(i, C));
        }
    }
    
    float d_val;

    while (!heap.empty()) {
        d_val = heap.top().d_i;
        cout << d_val << endl;
        heap.pop();
    }

    
    // main loop
    // invariant - while T^{*} is not of size n-1 
    
    //while (t_star->edges.size() < n-1) {
    //}
    
    return;
}

float min_ratio_st(int n, int m, vector<int>& neighbors, vector<int>& c, vector<int>& tau, vector<int>& indexes){
}

int main(int argc, char* argv[]) {
    const clock_t time_0 = clock();

    // const std::string datafile = argv[1];
    // const std::string outputfile = argv[2];
    
    int n = 3;
    int m = 3;

    // test graph G = (V, E)
    // V = {0, 1, 2}
    // E = {(0, 1), (0, 2), (1, 2)}
    // c = tau = 1

    vector<int> neighbors = {1, 2, 0, 2, 0, 1};
    vector<int> c = {1, 1, 1, 1, 1, 1};
    vector<int> tau = {1, 1, 1, 1, 1, 1};
    vector<int> indexes = {0, 2, 4};
    vector<float> h = {1, 1, 1, 1, 1, 1};
    float C = 1.0;
    MST t_star = MST();

    mst_prims(n, m, C, &neighbors, &h, &indexes, &t_star);

    // read graph file 
    // assign n,m
    // one continous list for each edge and costs, so 3 vectors total
    // neighbors (size 2m)
    // c (size 2m)
    // tau (size 2m)
    // use vector indexes (size n) to store index starts for neighbor of node i

    // float obj = min_ratio_st(n, m, &neighbors, &c, &tau, &indexes)

    float run_time = float(clock() - time_0) / CLOCKS_PER_SEC;

    return 0;
}
