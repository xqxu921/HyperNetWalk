"""
functions for permuting a mutation profile
"""
import random
import networkx as nx
import logging
# logging.basicConfig(level=logging.DEBUG)

from . import misc


def bipartite_double_edge_swap(G, genes, samples, nswap=1, max_tries=1e75):
    """A modified version of bipartite_double_edge_swap function from multi Dendrix,
    which is a modified version of double_edge_swap in NetworkX
    preserve the bipartite structure of the graph.

    :param G: nx.Graph B(G, S) a bipartite graph
    :param genes: list of genes
    :param samples: list of samples
    :param nswap: int, number of double edge swap to perform
    :param max_tries:int, maximum number of attests to swap edges
    :return: nx.Graph, permuted graph
    """
    if G.is_directed():
        raise nx.NetworkXError("double_edge_swap() not defined for directed graphs.")
    if nswap > max_tries:
        raise nx.NetworkXError("Number of swaps > number of tries allowed.")
    if len(G) < 4:
        raise nx.NetworkXError("Graph has less than four nodes.")

    # Instead of choosing uniformly at random from a generated edge list,
    # this algorithm chooses nonuniformly from the set of nodes with
    # probability weighted by degree.

    n=0
    swapcount=0

    gkeys,gdegrees=zip(*G.degree(genes)) # keys, degree for genes
    gcdf=nx.utils.cumulative_distribution(gdegrees)  # cdf of degree for genes

    pkeys,pdegrees=zip(*G.degree(samples)) # keys, degree for samples
    pcdf=nx.utils.cumulative_distribution(pdegrees)  # cdf of degree for samples

    while swapcount < nswap:
        # pick two random edges without creating edge list
        # choose source node indices from discrete distribution
        gi=nx.utils.discrete_sequence(1,cdistribution=gcdf)
        pi=nx.utils.discrete_sequence(1,cdistribution=pcdf)

        gene1=gkeys[gi[0]] # convert index to label
        sample1=pkeys[pi[0]]

        sample2 = random.choice(list(G[gene1]))
        gene2 = random.choice(list(G[sample1]))

        # don't create parallel edges
        if (gene1 not in G[sample1]) and (gene2 not in G[sample2]):
            G.add_edge(gene1,sample1)
            G.add_edge(gene2,sample2)

            G.remove_edge(gene1,sample2)
            G.remove_edge(gene2, sample1)
            swapcount+=1
        if n >= max_tries:
            e=('Maximum number of swap attempts (%s) exceeded '%n +
            'before desired swaps achieved (%s).'%nswap)
            raise nx.NetworkXAlgorithmError(e)
        n+=1
        if n % 10000 == 0:
            logging.debug("%d swaps..\n" %n)
    return G


def permute_mut_graph(G, genes, samples, Q=100):
    """Permutes a given mutation profile B(G, S) by performing |E| * Q edge swaps.

    :param G: nx.Graph B(G, S) a bipartite graph
    :param genes: list of genes
    :param samples: list of samples
    :param Q: constant multiplier for number Q * | E | of edge swaps to perform (default and suggested value: 100).
    See `Milo et al. (2003) <http://goo.gl/d723i>`_ for details on choosing Q.

    :returns: H: nx.Graph permuted bipartite graph
    """

    H = G.copy()
    H = bipartite_double_edge_swap(H, genes, samples, nswap=Q * len(G.edges()))
    return H


def construct_mut_graph_per_type(mut_dic, cancers, type_idx_dic):
    """ given mutation profile between genes and samples,
    create a bipartite graph for each cancer type separately

    :param mut_dic: mut dictionary gene -> list of cover weights for each sample
    :param cancers: cancer subtypes
    :param type_dic: dict cancer -> sample indices
    :return mut_graphs: alteration bipartite graphs for each cancer type
                dict cancer -> nx.Graph
    """
    mut_graphs = {}
    for cancer in cancers:
        mut_graphs[cancer] = nx.Graph()
    for gene in mut_dic:
        # mutated samples
        mut_sam_idxs = misc.get_positives(mut_dic[gene])
        for cancer in cancers:
            # create edges to mutated samples in a given cancer type
            edges = [(gene, s) for s in set(type_idx_dic[cancer]).intersection(mut_sam_idxs)]
            mut_graphs[cancer].add_edges_from(edges)
    return mut_graphs


def construct_mut_graph_per_type_from_list(mut_list_dic, cancers, type_idx_dic):
    """ given mutation profile "AS A LIST FORMAT" between genes and samples,
    create a bipartite graph for each cancer type separately

    :param mut_list_dic: mut dictionary gene -> list of covered sample idxs
    :param cancers: cancer subtypes
    :param type_dic: dict cancer -> sample indices
    :return mut_graphs: alteration bipartite graphs for each cancer type
                dict cancer -> nx.Graph
    """
    mut_graphs = {}
    for cancer in cancers:
        mut_graphs[cancer] = nx.Graph()
    for gene in mut_list_dic:
        # mutated samples
        mut_sam_idxs = mut_list_dic[gene]
        for cancer in cancers:
            # create edges to mutated samples in a given cancer type
            edges = [(gene, s) for s in set(type_idx_dic[cancer]).intersection(mut_sam_idxs)]
            mut_graphs[cancer].add_edges_from(edges)
    return mut_graphs


def construct_mut_dic_from_graphs(graphs, genes, nsamples):
    """  construct mutation dic from permuted bipartite graphs for all types

    :param 	graphs:  cancer type -> a bipartite graph B(G, S)
                S is given as sample indices to make it easy to construct the mutation list
    :param 	genes: list of genes
    :param 	nsamples: number of samples
    :return
            mut_dic: dict gene -> altered or not for each sample (in the order as in samples)
    """
    mut_dic = {}
    for gene in genes:
        mut_dic[gene] = [0 for i in range(nsamples)]
        for cancer in graphs:
            if gene not in graphs[cancer]:
                continue
            for i in graphs[cancer].neighbors(gene):
                mut_dic[gene][i] = 1
    return mut_dic
