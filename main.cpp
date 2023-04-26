#include <iostream>
#include "graph.h"
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>

void getDatasets(const string fname, vector<string>& names) {
    string currentDataset;
    std::ifstream datasetFiles(fname);
    while (datasetFiles >> currentDataset) {
        names.push_back(currentDataset); 
    }
    return;
}

int main() {
    vector<string> dataset_names;
    int s_nodes;
    double s_density;
    ofstream runtimes_file;
    runtimes_file.open("algo2_amazon_runtime2.txt");
    getDatasets("amazon_datasets/amazon_datasets.txt", dataset_names);
    for (const string s : dataset_names) {
        runtimes_file << s << "\n";
	string edgefile = "amazon_datasets/" + s + ".edges";
        string colorfile = "amazon_datasets/" + s + ".colors";
        auto [G, colors] = readData(edgefile, colorfile);
        double alpha_values[7] = {0.5,0.55,0.6,0.7,0.8,0.9,0.99};
        for (auto a : alpha_values) {
            auto start = std::chrono::high_resolution_clock::now();
            auto [s_nodes, s_density] = algorithm2(G, colors, a);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            runtimes_file << a << " " << s_density << " " << s_nodes.size() << " " << duration.count() << "\n";
        }
    }
    runtimes_file.close();
    return 0;
}
