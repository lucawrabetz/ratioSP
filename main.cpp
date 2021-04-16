#include <iostream>
#include <time.h>
#include <vector>
#include <fstream>
#include <queue>

using std::vector;

struct MST {
public:
    float obj_val; 
    vector<int> edges;

    MST() {
    }

    void add_edge(int _edge_index) {
        edges.push_back(_edge_index);
    }
    
    void set_objective(float _new_obj_val) {
        obj_val = _new_obj_val;
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
};


vector<int> mst_prims(int n, int m, vector<int>& neighbors,vector<int>& h, vector<int>& indexes) {
    /*
     * Subroutine for finding MST(k)
     * 
     * Receive current graph topology and costs 
     * Initialize --> heap, distance labels, pred, T^{*} (empty)
     * While T^{*} is not n-1, pop next node, update heap for adjacent 
     */
    
    // pre initialization
    // obj_val is total weight of tree
    // t_star is the MST (indexes of edges as per lists neighbors, h, c, tau)
    
        
    // initialize heap of vertices 
    // min_node has d(i) = 0 (we don't include in heap - first iteration pops and adds best neighbor of i)
    // neighbors of min_node have d(j) = h_{ij}
    // everyone else has d(j) = C,



    // main loop
    // invariant - while T^{*} is not of size n-1 
    
    while (t_star.size() < n-1) {
    }
    
    // return optimal tree (size n, last value is objective)
    t_star.push_back(obj);
    return t_star;
}

float min_ratio_st(int n, int m, vector<int>& neighbors, vector<int>& c, vector<int>& tau, vector<int>& indexes){
}

int main(int argc, char* argv[]) {
    const clock_t time_0 = clock();

    const std::string datafile = argv[1];
    const std::string outputfile = argv[2];
    
    int n;
    int m;

    vector<int> neighbors;
    vector<int> c;
    vector<int> tau;
    vector<int> indexes;

    // read graph file 
    // assign n,m
    // one continous list for each edge and costs, so 3 vectors total
    // neighbors (size m)
    // c (size m)
    // tau (size m)
    // use vector indexes (size n) to store index starts for neighbor of node i

    // float obj = min_ratio_st(n, m, &neighbors, &c, &tau, &indexes)

    float run_time = float(clock() - time_0) / CLOCKS_PER_SEC;

    return 0;
}
