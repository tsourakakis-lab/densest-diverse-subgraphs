#ifndef FAIRDENSESUBGRAPH_GRAPH_H
#define FAIRDENSESUBGRAPH_GRAPH_H

#include <queue>
#include <map>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <cmath>
using namespace std;

typedef vector<vector<int>> Graph;
typedef vector<int> Colormap;
typedef unordered_set<int> NodeSet;
typedef priority_queue<pair<double, int>, vector<pair<double, int>>, less<pair<double, int>>> maxHeap;

std::tuple<Graph,Colormap> readData(const string& edgelist, const string& colormap){
    Graph G;
    int numberNodes,numberEdges;
    int currentNodes = 0;
    int currentEdges = 0;
    unordered_map<int,int> relabeling;
    int src, dst;
    int first_line = 0;
    // Read graph from edge-file
    std::ifstream edgefile(edgelist);
    while (edgefile >> src >> dst) {
        if (first_line == 0) {
            numberNodes = src;
            numberEdges = dst;
            ++first_line;
            G.resize(numberNodes);
            continue;
        }
        if (relabeling.find(src) == relabeling.end()) {
            relabeling[src] = currentNodes++;
        }
        if (relabeling.find(dst) == relabeling.end()) {
            relabeling[dst] = currentNodes++;
        }
        G[relabeling[src]].push_back(relabeling[dst]);
        G[relabeling[dst]].push_back(relabeling[src]);
        currentEdges++;  
    }
    cout << currentNodes << "  " << numberNodes;
    assert(currentNodes == numberNodes);
    assert(currentEdges == numberEdges);

    // Read colors from file
    std::ifstream colorfile(colormap);
    Colormap nodes2colors;
    nodes2colors.resize(numberNodes);
    int node,color,nodeId;
    while(colorfile >> node >> color) {
        if (relabeling.find(node) == relabeling.end()) {
            // TODO: Raise error
        }
        nodeId = relabeling[node];
        nodes2colors[nodeId] = color;
    }
    printf("Read graph with %d nodes and %d edges\n", numberNodes, numberEdges);
    return {G,nodes2colors};
}

/* TODO: Re-write function to leverage unweighted graphs for better runtime
map and set
 */
NodeSet densestAtLeastK(const Graph& G, const int k) {
    int numberNodes = G.size();
    double numberEdges = 0.0; // Could be generalized degree
    auto *const degrees = new double[numberNodes];
    auto *const removalIndex = new int[numberNodes];
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> q;
    for (int i = 0; i < numberNodes; ++i) {
        degrees[i] = G[i].size();
        q.push(make_pair(degrees[i], i));
        numberEdges += degrees[i];
        removalIndex[i] = 0;
    }
    auto currentNodes = numberNodes;
    auto currentEdges = numberEdges / 2.0;
    auto bestSubgraphIdx = numberNodes;
    double currentDensity = currentEdges / (double) numberNodes;
    int curNode;
    while (!q.empty()) {
        curNode = q.top().second;
        q.pop();
        if (degrees[curNode] < 0) {
            continue;
        }
        removalIndex[curNode] = currentNodes;
        --currentNodes;
        currentEdges -= degrees[curNode];
        if ((currentEdges / (double) currentNodes) > currentDensity) {
            currentDensity = currentEdges / (double) currentNodes;
            bestSubgraphIdx = currentNodes;
        }
        // Update values for other nodes
        degrees[curNode] = -1;
        for (const auto neighbor : G[curNode]) {
            if (degrees[neighbor] > 0) { //TODO: Is there a way to save space here w/o using lazy deletion?
                q.push(make_pair(--degrees[neighbor], neighbor));
            }
        }
        if (numberNodes <= k) {
            break;
        }
    }
    NodeSet densestSubgraph;
    for (int i = 0; i < numberNodes; ++i) {
        if (removalIndex[i] <=  bestSubgraphIdx) {
            densestSubgraph.insert(i);
        }
    }
    std::cout << "Densest Subgraph Density: \n" << currentDensity << std::endl;
    return densestSubgraph;
}

///@brief : Function that returns the maximum element of a map
double getMaxFromVec(const vector<int>& v) {
    double maxElement = *std::max_element(v.begin(), v.end());
    return maxElement;
}

int getMinIdxVec(const vector<int>& v, const double maxColor, const vector<int>& outside_sizes, bool flag) {
    // Return idx of min as long as it is not equal to max and contains nodes
    int min_idx = -1;
    if (flag) {return min_idx;}
    double min_value = maxColor;
    for (int i=0; i < v.size(); ++i) {
        if ((v[i] < min_value) && (outside_sizes[i] > 0)) {
            min_idx = i;
            min_value = v[i];
        }
    }
    //if ((min_idx != -1) && (outside_sizes[min_idx] <= 0)) {min_idx=-1;} 
    return min_idx;
}

///@brief Returns worst candidate (one with less neighbors in subgraph)
int getWorstCandidate(const Graph& G, const NodeSet& candidates, const NodeSet& subgraph){
    int lessCn = G.size() + 1;
    int worstCandidate = -1;
    for (auto& c : candidates){
        int cn = 0;
        for (auto& neighbor : G[c]) {
            if (subgraph.find(neighbor) != subgraph.end()) {
                cn++;
            }
        }
        if (cn < lessCn) {
            lessCn = cn;
            worstCandidate = c;
        }
    }
    return worstCandidate;
}


// @brief Iteratively add nodes from color of minimum size, till constraints satisfied.
// If graph is not balanced, e.g. a(V)>a, then this may not work and we will also have
// to remove nodes.
void diversify(const Graph& G, NodeSet& S, const Colormap& colors, int numberOfColors, double a){
    // Keep nodes outside DS in max heaps, where (max) neighbors inside -> nodes_outside
    vector<int> color_sizes(numberOfColors,0); // Nodes with that color in S
    vector<int> outside_sizes(numberOfColors,0); // Nodes with that color not in S
    unordered_map<int,maxHeap> colors_outside; // list of nodes with that color not in S
    unordered_map<int,NodeSet> colors_inside; // list of nodes with that color not in S
    unordered_map<int,double> s_degree; // degree of nodes in S

    int numberOfNodes = G.size();
    for (int node=0; node < numberOfNodes; ++node){
        auto node_color = colors[node];
        if (S.find(node) != S.end()) {
            ++color_sizes[node_color];
            colors_inside[node_color].insert(node);
        } else {
            // Find neighbors in S and add size to max heap
            int nodes_in_s = 0;
            for (const auto& neighbor : G[node]) {
                if (S.find(neighbor) != S.end()) {
                    nodes_in_s += 1;
                }
            }
            s_degree[node] = nodes_in_s;
            colors_outside[node_color].push(make_pair(s_degree[node],node));
            ++outside_sizes[node_color];
        }
    }

    double maxColorSize;
    double sSize = S.size();
    bool pruneFlag = false;
    while (maxColorSize = getMaxFromVec(color_sizes), (maxColorSize / sSize) > a) {
        // as long as there are colors with size smaller than max
        // add nodes from those to subgraph
        auto min_color = getMinIdxVec(color_sizes, maxColorSize, outside_sizes, pruneFlag);
        if (min_color != -1) {
            auto best_node = colors_outside[min_color].top().second;
            colors_outside[min_color].pop();
            if (s_degree[best_node] < 0) continue;
            ++color_sizes[min_color];
            --outside_sizes[min_color];
            colors_inside[min_color].insert(best_node);
            S.insert(best_node);
            ++sSize;
            // Update heap for nodes outside S
            for (const auto neighbor : G[best_node]) {
                auto n_color = colors[neighbor];
                if (S.find(neighbor) == S.end()) {
                    colors_outside[n_color].push(make_pair(++s_degree[neighbor], neighbor));
                }
            }
            s_degree[best_node] = -1;
        } else {
            // Graph is not balanced, not possible to add nodes
            // We need to remove nodes from subgraph
            // Pick the ones from max color with less neighbors in S
            pruneFlag = true;
            auto max_color = std::max_element(color_sizes.begin(),color_sizes.end()) - color_sizes.begin();
            auto& candidate_nodes = colors_inside[max_color];
            auto worst_node = getWorstCandidate(G,candidate_nodes,S);
            S.erase(worst_node);
            candidate_nodes.erase(worst_node);
            --color_sizes[max_color];
            --sSize;
        }
    }

    return;
}

double getSubgraphDensity(const Graph& G, const NodeSet& S, const Colormap& colors) {
    auto numberOfNodes = S.size();
    double numberOfEdges = 0.0;
    unordered_map<int,int> colorCounts;
    for (const auto& u : S) {
        auto u_color = colors[u];
        ++colorCounts[u_color];
        for (const auto& v : G[u]) {
            if (S.find(v) != S.end()) {
                ++numberOfEdges;
            }
        }
    }
    numberOfEdges /= 2.0;
    double maxSize = -1;
    double totalSize = 0;
    for (auto const& c : colorCounts) {
        if (c.second > maxSize){
            maxSize = c.second;
        }
        totalSize += c.second;
    }
    cout << "Alpha: " << maxSize/totalSize << "\n";
    return numberOfEdges / (double) numberOfNodes;
}

int getNumberColors(const Colormap& colors) {
    unordered_set<int> unique_values;
    for (int c : colors) {
        unique_values.insert(c);
    }
    return unique_values.size();
}

std::tuple<NodeSet,double> algorithm2(const Graph& G, const Colormap& colors, const double a){
    auto k = ceil(1/ (double) a);
    auto S = densestAtLeastK(G,k);
    auto subgraph_density = getSubgraphDensity(G,S,colors);
    printf("Dalk Subgraph Density: %f Nodes: %ld\n", subgraph_density,S.size());
    int numberOfColors = getNumberColors(colors);
    cout << "Number of colors: " << numberOfColors << "\n";
    diversify(G,S,colors,numberOfColors,a);
    subgraph_density = getSubgraphDensity(G,S,colors);
    printf("After diversification density: %f Nodes: %ld \n", subgraph_density, S.size());
    return {S,subgraph_density};
}

///@brief Returns true if subgraph is connected and
/// false otherwise.
bool restrictedBFS(const Graph& G, const NodeSet& subgraph){
    queue<int> q;
    NodeSet visited;
    auto src = *subgraph.begin();
    q.push(src);
    visited.insert(src);
    while (!q.empty()) {
        src = q.front();
        q.pop();
        for (auto u : G[src]){
            if (subgraph.find(u) != subgraph.end()) {
                if (visited.find(u) == visited.end()) {
                    q.push(u);
                    visited.insert(u);
                }
            }
        }
    }
    return subgraph.size() == visited.size();
}

#endif //FAIRDENSESUBGRAPH_GRAPH_H
