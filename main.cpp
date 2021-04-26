#include <algorithm>
#include <iostream>
#include <time.h>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <queue>
#include <string>
#include <cmath>

using std::set;
using std::vector;
using std::priority_queue;
using std::cout;
using std::cin;
using std::endl;
using std::sort;
using std::floor;

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
    
    float getH() const {
        return H;
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

class EdgeComparator
{
public:
    float operator() (const Edge& e1, const Edge& e2)
    {
        return e1.getH() > e2.getH();
    }
};


//void mst_kruskal(int n, int m, float C, const vector<int>& neighbors, const vector<int>& indexes, const vector<int>& c, MST* t_star) {
//    /*
//     * Subroutine for finding MST(k) - kruskals
//     * 
//     * Receive current graph topology and costs 
//     * Receive MST address to update
//     * Initialize --> heap of edges and costs, and the bool vector of whos in the tree already (nodes)  
//     * While T^{*}.edges is not has size less than n-1, pop next edge, check for cycle, add edge 
//     * Updates the MST struct T^{*} without returning anything 
//     */
//    
//    // pre initialization
//    // t_star is the MST (indexes of edges as per lists neighbors, h, c, tau)
//    // t_star has attributes edges and obj_val
//    // other values - , 
//    t_star->clear(); 
//    
//        
//    // initialize heap of edges 
//    priority_queue <Edge, vector<Edge>, EdgeComparator> heap;
//    //vector<int> pred; // vector of predecessors
//    vector<bool> S; // who is already in the tree? (nodes)
//    int temp_index; // to hold the vertex index
//    float temp_dval; // to hold the distance label
//    int start_index; // to hold the index to start looping through in edge list
//    int stop_index; // to hold the index to stop looping through in edge list
//
//    int min_edge_index; // for the index of the minimum edge to be added to t_star
//    float min_edge_cost; // for the cost of the minimum edge to be added to t_star
//
//    for (int i = 0; i < n; ++i) {
//        S.push_back(false);
//    }       
//
//    for (int e = 0; e < m; ++e) {
//    }
//
//    start_index = indexes[0]; 
//    stop_index = indexes[1];
//
//    for (int j = start_index; j < stop_index; ++j) {
//        // push neighbors of node 0 to heap with distance label c_0j and set their pred = 0
//        heap.push(Vertex(neighbors[j], c[j])); // add node 0 to heap 
//        pred[neighbors[j]] = 0;
//    }    
//
//    // set S[0] to true - consider node 0 in the tree 
//    S[0] = true; 
//
//    /* use this loop to test that the heap is initializing correctly 
//    float d_val;
//    while (!heap.empty()) {
//        d_val = heap.top().d_i;
//        cout << d_val << endl;
//        heap.pop();
//    }
//    */
//
//    // main loop
//    // invariant - while T^{*} is not of size n-1 
//    
//    while (t_star->edges.size() < n-1) {
//        // get the vertex with next lowest d_{i} from the heap, and delete from the heap
//        temp_index = heap.top().i; 
//        temp_dval = heap.top().d_i;
//        heap.pop();
//
//        start_index = indexes[temp_index];
//        stop_index = indexes[temp_index+1];
//
//        for (int e = start_index; e < stop_index; ++e) {
//            if (!S[neighbors[e]]) {
//                // if the neighbor is not already in S, push neighbor of node to heap with distance label c and set their pred
//                heap.push(Vertex(neighbors[e], c[e])); // add node 0 to heap 
//                pred[neighbors[e]] = temp_index;
//            }
//            if (neighbors[e] == pred[temp_index]) {
//                // set the index of the edge we are adding to e - this way we can keep track of it t_star by edge index and not as a tuple
//                min_edge_index = e;
//                min_edge_cost = c[e];
//            }
//        }
//        
//        t_star->add_edge(min_edge_index);
//        t_star->set_objective(t_star->obj_val + min_edge_cost);
//        // update the heap d_i values of each neighbor - !!!! we will re-add these neighbors 
//        // since priority_queue has no decrease_key implementation, so we will have duplicates
//        // but the ones with d = C+1 will not get popped so its fine
//        
//    }
//    
//    return;
//}


void mst_prims(int n, int m, float C, const vector<int>& neighbors, const vector<int>& indexes, const vector<float>& c, MST* t_star) {
    /*
     * Subroutine for finding MST(k) - prims
     * 
     * Receive current graph topology and costs 
     * Receive MST address to update
     * Initialize --> heap of nodes and distance labels, pred, clear T^{*} 
     * While T^{*}.edges is not has size less than n-1, pop next node, update heap for adjacent 
     * Updates T^{*}, doesn't return anything
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
    vector<float> best_ds; // keep basically a copy of the heap with the current 
    // (continued from above) distance label for every node
    // we need it because we don't decrease_key, we just add copies, 
    // and we need to be able to check the current values of things in the heap for the algo
    int temp_index; // to hold the vertex index
    float temp_dval; // to hold the distance label
    int start_index; // to hold the index to start looping through in edge list
    int stop_index; // to hold the index to stop looping through in edge list

    int min_edge_index; // for the index of the minimum edge to be added to t_star
    float min_edge_cost; // for the cost of the minimum edge to be added to t_star

    for (int i = 0; i < n; ++i) {
        pred.push_back(-1);
        S.push_back(false);
        best_ds.push_back(C+1);
    }       

    start_index = indexes[0]; 
    stop_index = indexes[1];

    for (int j = start_index; j < stop_index; ++j) {
        // push neighbors of node 0 to heap with distance label c_0j and set their pred = 0
        heap.push(Vertex(neighbors[j], c[j])); // add node 0 to heap 
        best_ds[neighbors[j]] = c[j];
        pred[neighbors[j]] = 0;
    }    

    // set S[0] to true - consider node 0 in the tree 
    S[0] = true; 

    //use this loop to test that the heap is initializing correctly 
    //float d_val;
    //while (!heap.empty()) {
    //    d_val = heap.top().d_i;
    //    cout << d_val << endl;
    //    heap.pop();
    //}
    

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
                heap.push(Vertex(neighbors[e], c[e])); // add node to heap 
                if (best_ds[neighbors[e]]>c[e]) {
                    // we only change the pred for this guy if this edge is his best edge! he may have a better edge from a previously set predecessor
                    pred[neighbors[e]] = temp_index;
                    best_ds[neighbors[e]] = c[e];
                }
            }
            if (neighbors[e] == pred[temp_index]) {
                // set the index of the edge we are adding to e - this way we can keep track of it t_star by edge index and not as a tuple
                min_edge_index = e;
                min_edge_cost = c[e];
            }
        }
        
        t_star->add_edge(min_edge_index);
        t_star->set_objective(t_star->obj_val + min_edge_cost);
        S[temp_index] = true;
        // update the heap d_i values of each neighbor - !!!! we will re-add these neighbors 
        // since priority_queue has no decrease_key implementation, so we will have duplicates
        // but the ones with d = C+1 will not get popped so its fine
    }
    
    return;
}

vector<float> get_k_ratios(int m, const vector<int>& c, const vector<int>& tau){
    // compute and sort ratios ce-cf/de-df for each pair of edges
    // receive vector of edge costs (indexed e = 1-m) c and tau
    //
    vector<float> k_ratios;
    float k_ratio;
    // TODO : check for duplicates and correct index at insertion
    for (int e = 0; e < m-1; ++e) {
        for (int f = e+1; f < m; ++f) {
            if (tau[e]-tau[f] == 0) {
                continue;
            }
            k_ratio = (c[e]-c[f])/(tau[e]-tau[f]);
            k_ratios.push_back(k_ratio);
        }
    }
    set<int> s;
    unsigned size = k_ratios.size();
    for( unsigned i = 0; i < size; ++i ) s.insert( k_ratios[i] );
    k_ratios.assign( s.begin(), s.end() );
    sort(k_ratios.begin(), k_ratios.end()); 
    return k_ratios;
}

float eval_tree(float k_ratio, const vector<int>& c, const vector<int>& tau, MST* tree) {
    /*
     * Return the objective value of the tree object on the parametrized H(k) problem
     */

    vector<int> edges = tree->edges;
    float h_curr;
    float h_tot = 0;
    for (int i = 0; i < edges.size(); ++i) {
        h_curr = c[edges[i]] - k_ratio*tau[edges[i]];
        h_tot += h_curr;
    }
    return h_tot;
}

float eval_tree_nonparam(const vector<int>& c, const vector<int>& tau, MST* tree) {
    /*
     * Return the objective value of the tree object on the parametrized H(k) problem
     */
    float obj; 
    int c_tot = 0;
    int tau_tot = 0;
    vector<int> edges = tree->edges;
    for (int i = 0; i < edges.size(); ++i) {
        c_tot += c[edges[i]];
        tau_tot += tau[edges[i]];
    }
    obj = c_tot / tau_tot;
    return obj;
}

float min_ratio_st(int n, int m, float C, const vector<int>& c_raw, const vector<int>& tau_raw, const vector<int>& raw_indexes, const vector<int>& neighbors, const vector<int>& c, const vector<int>& tau, const vector<int>& indexes, MST* tree){
    /*
     * Subroutine for finding MRST 
     * Receive current graph topology and costs 
     * Initialize -->  
     * Anything with a _t is updated every iteration
     */
    
    // initialize 
    int t = 0;
    float k_t;
    float obj_val = 0;
    
    float A_val;
    float B_val;
    vector<float> k_vals = get_k_ratios(m, c_raw, tau_raw); 
    int r = k_vals.size();
    int alpha_t = 0;
    int beta_t = r;
    int j_index = 0;
    vector<float> h;
    h.reserve(c.size());

    // log k_ratios
    // cout<< "k ratios: " << endl;
    // for (int i = 0; i < r; ++i){
    //     cout<< k_vals[i] << endl;
    // }

    bool optimal = false;

    while (!optimal) {
        // set j index and k_t 
        j_index = floor((alpha_t+beta_t)/2);
        k_t = (k_vals[j_index]+k_vals[j_index+1])/2;

        // set the h_vector
        for (int e = 0; e < c.size(); ++e) {
            h[e] = c[e] - k_t * tau[e];
        }

        // resolve the tree object on this h function
        mst_prims(n, m, C, neighbors, indexes, h, tree);

        // calculate A_val and B_val to update everything

        optimal = true;
    }

    return obj_val;
}

int main(int argc, char* argv[]) {
    // cout << "main started" << endl; 
    const clock_t time_0 = clock();
    
    // to use the following test graph G = (V, E) pass 'complete3.graph' as datafile
    // V = {0, 1, 2}
    // E = {(0, 1), (0, 2), (1, 2)}
    // c = tau = 1
    // when you have args uncomment and use the next two lines const std::string datafile = argv[1]; const std::string outputfile = argv[2]; to hard code the datafiles const std::string datafile = "example.graph";
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
    vector<int> c_raw;
    vector<int> tau_raw;
    vector<int> c;
    vector<int> tau;
    vector<int> indexes;
    vector<int> raw_indexes;
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
                    c_raw.push_back(edge_c);
                    tau_raw.push_back(edge_tau);
                    raw_indexes.push_back(counter-2);
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
            cout<<j<<": ("<<std::to_string(neighbors[j]);
            cout<<","<<std::to_string(c[j]);
            cout<<","<<std::to_string(tau[j])<<")";
            cout<<endl;
        }
    }


    const vector<int> final_neighbors = neighbors;
    const vector<int> final_c = c;
    const vector<int> final_tau = tau;
    const vector<int> final_c_raw = c_raw;
    const vector<int> final_tau_raw = tau_raw;
    const vector<int> final_raw_indexes = raw_indexes;
    const vector<int> final_indexes = indexes;
    MST* t_star = new MST();
    

    float obj_val = min_ratio_st(n, m, final_c_raw, final_tau_raw, final_raw_indexes, final_neighbors, final_c,  final_tau, final_indexes);
    // mst_prims(n, m, C, final_neighbors, final_indexes, final_c, t_star);
    //cout<<"prims obj: "<<t_star->obj_val<<endl;
    //for (int i = 0; i < n-1; ++i) {
    //    cout<<t_star->edges[i]<<endl;
    //}



    float run_time = float(clock() - time_0) / CLOCKS_PER_SEC;

    return 0;
}
