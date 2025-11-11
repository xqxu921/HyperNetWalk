"""
various io utility functions for processing TCGA data
read/write etc.
"""

import sys
import pandas
import mut_ex.misc


def read_dic(filename):
    """ read a file with two columns and create a dic
        the first column as keys and the second as values
    :param filname
    :return: dict first column -> the second column
    """
    lines = open(filename).readlines()
    data_dic = {}
    for l in lines:
        tkns = l.split()
        if len(tkns) > 2:
            sys.stderr.write("more than two columns")
            return {}
        data_dic[tkns[0]] = tkns[1]
    return data_dic


def read_mut_list(filename, sep=",", samples=False, genes=True):
    """
    read a bipartite graph B(G, S) in the form of list of mutated samples
        gene \t cs where cs for comma separated indices (default)
        G: genes, S: samples
    ** USE THIS COMPACT FORMAT ONLY FOR UNWEIGHTED CASES.
    example:
        PTEN	12,14,36,40
    :param filename
    :param sep default=","
    :param samples True if the first row is list of samplaes
    :param genes True if the first column is gene name
    :return
            if genes is True:
                dict gene -> list of sample indices
                if samples is True:
                    slist --> list of sample names
            elif genes is False:
                list of sample indices list
                if samples is True:
                    slist --> list of sample names
    """
    # todo: check output format usage change:  sample list is returned separately
    # create gene->list
    if genes is True:
        temp_dic = read_dic(filename)
        rel_dic_list = {}
        if samples is True:
            slist = temp_dic["samples"].split(sep) # samples
            del temp_dic["samples"]
        else:
            slist = []
        for g in temp_dic:
            if len(temp_dic[g]) == 0: # remove the entry if the gene is not mutated in any samples
                continue
            rel_dic_list[g] = [int(x) for x in temp_dic[g].split(sep)]
        return rel_dic_list, slist
    # create list of lists
    elif genes is False:
        idx_list = []
        lines = open(filename).readlines()
        if samples is True:
            slist = lines[0].strip().split(sep)
        else:
            idx_list.append([int(x) for x in lines[0].strip().split(sep)])
            slist = []
        for l in lines[1:]:
            idx_list.append([int(x) for x in l.strip().split(sep)])
        return idx_list, slist


def read_mut_matrix(filename, mw=3, mutsig_file=None):
    """ read a bipartite graph B(G, S) in the form of a labeled matrix
        G: genes, S: samples
        the first row is the list of sample labels
        the first column is for gene names
        Each element is labeled as N(None), A(Amp), D(Del), M(Somatic Mutation), or AM/DM
        example:
                    s1	s2	s3	s4 	S5  ...
            PTEN	N 	C	N	M	B   ...

        labels are converted to weights based on weight_dic (if mw > 0)
        weight_dic = dict([('N', 0), ('A', 1), ('D', 1), ('M', mw), ('AM', mw+1), ('DM', mw+1)])
        if mutsig_file is given, multiply mutsig score + 1 to the each scores (IGNORE mw)

    :param filename
    :param mw (the relative weight of somatic mutation) default=3
            if mw < 0 , no conversion
    :param mutsig_file: gene weight file from mutsig
        ex.
            1 gene mutsig_score
            2 ABCA1 -0.0
            3 ABCA2 0.0862616260817

    :return genes: list of genes
        samples: list of samples
        data_dic: dict gene g -> the list (possibly) converted from the label to edge weight e(g, s)
    """
    lines = open(filename).readlines()
    genes = []
    samples = lines[0].split()[1:]
    data_dic = {}

    if mutsig_file is not None:  # if mutsig file is given
        mutsig = pandas.read_table(mutsig_file, sep=" ", index_col=0).to_dict()['mutsig_score']
        mutsig_data_dic = {}
        for g in data_dic:
            # multiply mutsig_score to the gene
            mutsig_data_dic[g] = [(mutsig[g]+1)*x for x in data_dic[g]]
        data_dic = mutsig_data_dic

    elif mw > 0:  # convert based on weight_dic/mw
        # weight_dic = dict([('N', 0), ('C', 1), ('M', mw), ('B', mw+1)])
        weight_dic = dict([('N', 0), ('A', 1), ('D', 1), ('M', mw), ('AM', mw+1), ('DM', mw+1)])
        for l in lines[1:]:
            tkns = l.split()
            data_dic[tkns[0]] = [float(weight_dic[x]) for x in tkns[1:]]
            genes.append(tkns[0])
    else:  # no conversion
        for l in lines[1:]:
            tkns = l.split()
            data_dic[tkns[0]] = [x for x in tkns[1:]]
            genes.append(tkns[0])

    return genes, samples, data_dic


def write_mut_matrix(genes, samples, data_dic, filename):
    """
    write a bipartite graph B(G, S) in the form of a weighted matrix
        G: genes, S: samples
        the first row is the list of sample names
        the first column is for gene names
        example:
                sample1\tsample2\t ...
            PTEN\t0\t1\t\t0\t3..

    :param genes: list of genes
    :param samples: list of samples
    :param data_dic: dict gene -> weights of edges in B(G, S), 0 for no edge
    :param filename
    :return
    """
    f = open(filename, 'w')
    f.write("gene\t%s\n" % "\t".join(samples))
    for g in genes:
        f.write("%s\t" %g)
        f.write("%s\n" % "\t".join([str(x) for x in data_dic[g]]))
    f.close()


def write_mut_list(rel_dic, filename, samples=None, sep=","):
    """
    write a bipartite graph B(G, S) in the form of list of positive samples
        gene\tcs where csi for comma separated indices (default)
        G: genes, S: samples
    ** USE THIS FORMAT FOR UNWEIGHTED CASES. COMPACT FORMAT
    example:
        PTEN	12,14,36,40
    :param rel_dic: dict gene -> weights of edges in B(G, S), 0 for no edge
    :param filename
    :param sep default=","
    :return
    """

    list_dic = {}
    for g in rel_dic:
        covers = mut_ex.misc.get_positives(rel_dic[g])
        if len(covers) == 0:
            continue
        list_dic[g] = sep.join([str(x) for x in covers])

    if samples is None:
        write_dic(list_dic, filename)
    else:
        write_dic(list_dic, filename, "samples\t"+",".join(samples))


def write_alt_list(list_dic, filename, samples=None, sep=","):
    """
    write list of sample indices
        in the format of
        gene\tcs where csi for comma separated indices (default)
        G: genes, S: samples
    ** USE THIS FORMAT FOR UNWEIGHTED CASES. COMPACT FORMAT
    example:
        PTEN	12,14,36,40
    :param list_dic: gene -> list of altered sample indices
    :param filename:
    :param samples: if sample is given, write in the first row
    :param sep: separator
    :return
    """
    for g in list_dic:
        list_dic[g] = sep.join([str(x) for x in list_dic[g]])

    if samples is None:
        write_dic(list_dic, filename)
    else:
        write_dic(list_dic, filename, "samples\t"+",".join(samples))


def write_dic(any_dic, filename, labels=None):
    """ write a file with two columns and create a dic
        the first column as keys and the second as values
    :param any_dic
    :param filename
    :return
    """

    f = open(filename, 'w')
    if labels is not None:
        f.write("%s\n" %labels)
    for x in any_dic:
        f.write("%s\t%s\n" % (str(x), str(any_dic[x])))
    f.close()


def create_type_idx(samples, types, sample_type_dic):
    """
    create a dict mapping type -> list of sample indices for the type
    the indices given in the samples are used

    :param samples: list of samples
    :param types: list of cancer types
    :param sample_type_dic: sample->type dic
    :return:
    """
    type_idx_dic = {}
    for ty in types:
        type_idx_dic[ty] = filter(lambda i: sample_type_dic[samples[i]] == ty, range(len(samples)))

    return type_idx_dic


