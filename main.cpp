#include <iostream>
#include <time.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <queue>
#include <string>

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
    vector<int> neighbors;
    vector<int> c;
    vector<int> tau;

    Vertex(int _i, float _d_i) {
        i = _i;
        d_i = _d_i;
    }

    void update_Vertex(float _new_d_i) {
        d_i = _new_d_i;
    }

    void update_neighbors(int _new_neighbor, int _new_c, int _new_tau) {
        neighbors.push_back(_new_neighbor);
        c.push_back(_new_c);
        tau.push_back(_new_tau);
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


void mst_prims(int n, int m, float C, const vector<int>& neighbors, const vector<int>& indexes, const vector<int>& c, MST* t_star) {
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
    vector<int> pred; // vector of predecessors
    vector<bool> S; // who is already in the tree? (nodes)
    int temp_index; // to hold the vertex index
    float temp_dval; // to hold the distance label
    int start_index; // to hold the index to start looping through in edge list
    int stop_index; // to hold the index to stop looping through in edge list

    int min_edge_index; // for the index of the minimum edge to be added to t_star
    float min_edge_cost; // for the cost of the minimum edge to be added to t_star

    for (int i = 1; i < n; ++i) {
        pred.push_back(-1);
        S.push_back(false);
    }       

    start_index = indexes[0]; 
    stop_index = indexes[1];

    for (int j = start_index; j < stop_index; ++j) {
        // push neighbors of node 0 to heap with distance label c_0j and set their pred = 0
        heap.push(Vertex(neighbors[j], c[j])); // add node 0 to heap 
        pred[neighbors[j]] = 0;
    }    

    // set S[0] to true - consider node 0 in the tree 
    S[0] = true; 

    /* use this loop to test that the heap is initializing correctly 
    float d_val;
    while (!heap.empty()) {
        d_val = heap.top().d_i;
        cout << d_val << endl;
        heap.pop();
    }
    */

    // main loop
    // invariant - while T^{*} is not of size n-1 
    
    while (t_star->edges.size() < n-1) {
        // get the vertex with next lowest d_{i} from the heap, and delete from the heap
        temp_index = heap.top().i; 
        temp_dval = heap.top().d_i;
        heap.pop();

        start_index = indexes[temp_index];
        stop_index = indexes[temp_index+1];

        for (int e = start_index; e < stop_index; ++e) {
            if (!S[neighbors[e]]) {
                // if the neighbor is not already in S, push neighbor of node to heap with distance label c and set their pred
                heap.push(Vertex(neighbors[e], c[e])); // add node 0 to heap 
                pred[neighbors[e]] = temp_index;
            }
            if (neighbors[e] == pred[temp_index]) {
                // set the index of the edge we are adding to e - this way we can keep track of it t_star by edge index and not as a tuple
                min_edge_index = e;
                min_edge_cost = c[e];
            }
        }
        
        t_star->add_edge(min_edge_index);
        t_star->set_objective(t_star->obj_val + min_edge_cost);
        // update the heap d_i values of each neighbor - !!!! we will re-add these neighbors 
        // since priority_queue has no decrease_key implementation, so we will have duplicates
        // but the ones with d = C+1 will not get popped so its fine
        
    }
    
    return;
}

float min_ratio_st(int n, int m, vector<int>& neighbors, vector<int>& c, vector<int>& tau, vector<int>& indexes){
}

int main(int argc, char* argv[]) {
    // cout << "main started" << endl; 
    const clock_t time_0 = clock();
    
    // to use the following test graph G = (V, E) pass 'complete3.graph' as datafile
    // V = {0, 1, 2}
    // E = {(0, 1), (0, 2), (1, 2)}
    // c = tau = 1
    // when you have args uncomment and use the next two lines
    // const std::string datafile = argv[1];
    // const std::string outputfile = argv[2];
    
    // to hard code the datafiles
    const std::string datafile = "complete3.graph";
    const std::string outputfile = "output.txt";
    
    // initialize file objects
    std::ifstream ifile (datafile);
    std::ofstream ofile (outputfile);
    vector<Vertex> vertices; // maintain a vector of vertices at all times - first to construct contiguous neighbor list, then for reference in heap creation
    // cout << "vertex vector created" << endl; 

    // read graph file 
    // assign n,m
    // one continous list for each edge and costs, so 3 vectors total
    // neighbors (size 2m)
    // c (size 2m)
    // tau (size 2m)
    // use vector indexes (size n) to store index starts for neighbor of node i

    int n;
    int m;

    vector<int> neighbors;
    vector<int> c;
    vector<int> tau;
    vector<int> indexes;
    vector<float> h;
    float C; 
    std::string line;
    
    // cout<<"all other vectors and info created"<<endl;

    if (ifile.is_open()) {
        int counter = 1;
        const char* cline;
        int u;
        int v;
        int edge_c;
        int edge_tau; 
        float max_c = 0;
        // cout<<"starting while loop on file"<<endl;
        while (getline(ifile, line)) {
            if (line != "") {
                if (counter > 1) {
                    // not first line - handle a new edge
                    cline = line.c_str();
                    sscanf(cline, "%d %d %d %d", &u, &v, &edge_c, &edge_tau);
                    vertices[u].update_neighbors(v, edge_c, edge_tau);
                    vertices[v].update_neighbors(u, edge_c, edge_tau);
                    if (edge_c > max_c) {
                        max_c = edge_c;
                    }
                }
                else {
                    // first line - assign n and m
                    cline = line.c_str();
                    sscanf(cline, "%d %d", &n, &m);
                    // create all vertex objects
                    for (int i = 0; i < n; ++i){
                        vertices.push_back(Vertex(i, 0));
                    }
                }
                ++counter;
            }
        }
        C = max_c;
    }
    
    // create contiguous list of edges
    int next_index_start = 0; 
    for (int i = 0; i<n; ++i) {
        indexes.push_back(next_index_start);
        for (int j = 0; j < vertices[i].neighbors.size(); ++j) {
            neighbors.push_back(vertices[i].neighbors[j]);
            h.push_back(vertices[i].c[j]); // temporarily use c for h for testing of MST subroutines
            c.push_back(vertices[i].c[j]);
            tau.push_back(vertices[i].tau[j]);
            ++next_index_start;
        }
    }
    // add an n+1th value to indexes - so the nth vertex has a 'stop value'
    indexes.push_back(2*m);
    
    
    /* uncomment to print graph after being read (by vertex objects)
    // for (int i = 0; i < n; ++i) {
    //     cout<<std::to_string(vertices[i].i)<<": neighbors (index, c, tau) - "<<endl; 
    //     for (int j = 0; j < vertices[i].neighbors.size(); ++j) {
    //         cout<<"("<<std::to_string(vertices[i].neighbors[j]);
    //         cout<<","<<std::to_string(vertices[i].c[j]);
    //         cout<<","<<std::to_string(vertices[i].tau[j])<<")";
    //         cout<<endl;
    //     }
    // }
    */

    // uncomment to print graph after being read (by full contiguous list)
    int start_index;
    int stop_index;
    for (int i = 0; i < n; ++i) {
        cout<<std::to_string(vertices[i].i)<<": neighbors (index, c, tau) - "<<endl; 
        start_index = indexes[i];
        stop_index = indexes[i+1];
        cout << "start index: " << std::to_string(start_index) << ", stop index: " << std::to_string(stop_index)<<endl;
        for (int j = start_index; j < stop_index; ++j) {
            cout<<"("<<std::to_string(neighbors[j]);
            cout<<","<<std::to_string(c[j]);
            cout<<","<<std::to_string(tau[j])<<")";
            cout<<endl;
        }
    }

    const vector<int> final_neighbors = neighbors;
    const vector<int> final_c = c;
    const vector<int> final_tau = tau;
    const vector<int> final_indexes = indexes;
    MST* t_star = new MST();
    
    mst_prims(n, m, C, final_neighbors, final_indexes, final_c, t_star);
    
    cout<<"prims obj: "<<t_star->obj_val<<endl;
    for (int i = 0; i < n-1; ++i) {
        cout<<t_star->edges[i]<<endl;
    }

    // float obj = min_ratio_st(n, m, &neighbors, &c, &tau, &indexes)

    float run_time = float(clock() - time_0) / CLOCKS_PER_SEC;

    return 0;
}
