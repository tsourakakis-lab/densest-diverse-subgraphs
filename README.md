# Densest diverse subgraphs: How to Plan a Successful Cocktail Party with Diversity
### Atsushi Miyauchi, Tianyi Chen, Kostas Sotiropoulos, Charalampos Tsourakakis 

## Description

Dense subgraph discovery methods are routinely used in a variety of applications including the identification of a team of skilled individuals for collaboration from a social network.  However, when the network's node set is associated with a sensitive attribute such as race, gender, religion, or political opinion, the lack of diversity can lead to lawsuits.  More broadly, the problem of finding diverse clusters appears in other domains including connectomics where neurons in the brain are associated with different  functional roles. 

In this work, we focus on the problem of finding densest diverse subgraphs in graphs whose nodes have different attribute values/types that we refer to as colors. We propose two novel formulations motivated by different realistic scenarios.    Our first formulation, called the  **densest diverse subgraph problem**, guarantees that no color represents more than an $\alpha$ fraction of the nodes in the output subgraph. By varying $\alpha$ we can  range the diversity constraint, and interpolate from diverse dense subgraphs  where all colors have to be equally represented to unconstrained dense subgraphs.  We design an $\Omega(1/\sqrt{n})$-approximation algorithm, where $n$ is the number of nodes, that is independent of the number of colors and scales gracefully to large-scale networks. The \ddsp generalizes and significantly improves in various ways the state-of-the-art due to Anagnostopoulos et al.. Our second formulation is motivated by the setting where the minority group is very small to the node set. We propose the {\it Dal$\vec{k}$S} problem, a novel generalization of the classic densest at least-$k$ subgraph problem (Dal$k$S), where instead of a single value $k$, we have a vector of cardinality demands $\vec{k}$ with one coordinate per color class. We design a $1/3$-approximation algorithm using linear programming together with an acceleration technique. We evaluate our  algorithms on synthetic and real-world datasets, including social and collaboration networks. Our proposed algorithms succeed in extracting  diverse clusters from social networks with sensitive node attributes.
