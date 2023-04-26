# Densest diverse subgraphs: How to Plan a Successful Cocktail Party with Diversity
### Atsushi Miyauchi, Tianyi Chen, Kostas Sotiropoulos, Charalampos E. Tsourakakis 

## Description

Dense subgraph discovery methods are routinely used in a variety of applications including the identification of a team of skilled individuals for collaboration from a social network.  However, when the network's node set is associated with a sensitive attribute such as race, gender, religion, or political opinion, the lack of diversity can lead to lawsuits.  More broadly, the problem of finding diverse clusters appears in other domains including connectomics where neurons in the brain are associated with different  functional roles. 

In this work, we focus on the problem of finding densest diverse subgraphs in graphs whose nodes have different attribute values/types that we refer to as colors. We propose two novel formulations motivated by different realistic scenarios.    Our first formulation, called the  **densest diverse subgraph problem (DDSP)**, guarantees that no color represents more than an $\alpha$ fraction of the nodes in the output subgraph. By varying $\alpha$ we can  range the diversity constraint, and interpolate from diverse dense subgraphs  where all colors have to be equally represented to unconstrained dense subgraphs.  We design an $\Omega(1/\sqrt{n})$-approximation algorithm, where $n$ is the number of nodes, that is independent of the number of colors and scales gracefully to large-scale networks. The **DDSP** generalizes and significantly improves in various ways the state-of-the-art due to Anagnostopoulos et al.. Our second formulation is motivated by the setting where the minority group is very small to the node set. We propose the DalkS problem, a novel generalization of the classic densest at least-k subgraph problem (DalkS), where instead of a single value k, we have a vector of cardinality demands **k** with one coordinate per color class. We design a $1/3$-approximation algorithm using linear programming together with an acceleration technique. We evaluate our  algorithms on synthetic and real-world datasets, including social and collaboration networks. Our proposed algorithms succeed in extracting  diverse clusters from social networks with sensitive node attributes.

## Instruction

### Dataset

Due to the Github file size limit, we only compress and upload partial Amazon networks in this Repository as "amazon_dataset.zip". Full datasets can be found from [1].

### Algorithm 2

- The code for Algorithm 2 is written in C++. All functions are included in the "graph.h" file.

* The function "readData(edgefile, colorfile)" takes as arguments the locations of two files:
** edgeFile: The input graph. The first line in the file is the number of vertices and the 
  number of edges. The remaining lines contain all edges in the graph. Notice that as the input
  graph is undirected all edges should be included just once.
** colorfile: The colors of the nodes in the graph. Each line is a pair of:
  "nodeId colorId", space-separated.
It returns a graph object (G) and a mapping <nodeId,colorId> (colors).

* The function "algorithm2(G, colors, a)" takes as its arguments the graph object (G) and the
  <nodeId,colorId> mapping of the previous function, as well the fairness parameter "a" and returns
  a set S of the nodes in the subgraph, and the density of it.

### Algorithm 6

- The Linear/Integer Programming of algorithm 6 is solved using Python scipy with highs solver. All functions as well as results are included in the "Alg6_exp" notebook.

[1] Aris Anagnostopoulos, Luca Becchetti, Adriano Fazzone, Cristina Menghini, and Chris Schwiegelshohn. 2020. Spectral Relaxations and Fair Densest Subgraphs. In Proceedings of the 29th ACM International Conference on Information & Knowledge Management (CIKM '20). Association for Computing Machinery, New York, NY, USA, 35–44. DOI:https://doi.org/10.1145/3340531.3412036
