"""
compute mutual exclusivity rank/pv via weighted sampling
"""

import numpy as np
import bisect

from . import misc


def compute_ex_cover(covers):
    """
    compute exclusive set
    :param covers: list of covers for a set of genes
    :return: ex_cover
    """
    ex_cov = set([])
    union_cov = set([])
    cross_cov = set([])
    for cov in covers:
        ex_cov = ex_cov.symmetric_difference(set(cov).difference(cross_cov))
        cross_cov = cross_cov.union(union_cov.intersection(cov))
        union_cov = union_cov.union(cov)
    return ex_cov


def compute_rank(cover_size, ws_cover_sizes, me_co="me", ordered=False):
    """
    given random instances (weighted sampling) and the original cover size,
    compute empirical pvalue for mutual exclusivity/co-occurence

    :param cover_size: original cover size
    :param ws_cover_sizes: list of random cover sizes
    :param me_co: me or co
    :param ordered: ws_cover_sizes is sorted if True
    :return rank: number of instances bigger than the original
    """

    if me_co == "me":
        if ordered:
            less = bisect.bisect_left(ws_cover_sizes, cover_size)
        else:
            less = len(list(filter(lambda c: c < cover_size, ws_cover_sizes)))
        rank = len(ws_cover_sizes)-less
    elif me_co == "co":
        if ordered:
            lesseq = bisect.bisect_right(ws_cover_sizes, cover_size)
        else:
            lesseq = len(list(filter(lambda c: c <= cover_size, ws_cover_sizes)))
        rank = lesseq
    return rank


def choose_random_tuples(ws_num, tsize, tnum, slack=0.3):
    """
    randomly choose patient indices tuples excluding redundant indices
    :param ws_num: num of random instances for each k in weighted sampling
    :param tsize: size of a tuple (number of genes)
    :param tnum: number of tuples to be chosen
    :param slack: fraction for overselection to avoid additional round of selection
    :return: rtuples: random indices of a tuple
    """
    rtuples = []
    cur_num = int(tnum*(1+slack))
    while True:
        ris = [np.random.choice(range(ws_num), cur_num) for i in range(tsize)]  # choose random indices
        cur_tuples = zip(*ris)
        rtuples.extend(filter(lambda ri: len(ri) == len(set(ri)), cur_tuples))  # no redundant indices in each tuple
        rtuples = list(set(rtuples))  # no redundant tuples
        if len(rtuples) >= tnum:
            return rtuples[:tnum]
        else: # not enough. choose more
            cur_num = int((tnum - len(rtuples))*(1+slack))
    return rtuples


def compute_ws_cover_sizes(cover_sizes, ws_k_cover_dic, tuple_num, ws_ex_cover_sizes_dic):
    """
    compute NULL cover sizes using weighted sampling
    check ws_ex_cover_sizes_dic first (reuse previous results if available)
    if not available or the number of samples not enough,
    sample more from ws_k_cover_dic and compute ws_cover_sizes
    and update ws_ex_cover_sizes_dic
    :param cover_sizes: cover sizes for a list of genes
    :param ws_k_cover_dic: k -> weighted sampling
    :param tuple_num: number of random ex_covers to compute
    :param ws_ex_cover_sizes_dic: to reuse the previous sampling
                (k tuples) -> list of ex cover sizes
    :return: ws_ex_cover_sizes: NULL ex cover sizes
             ws_ex_cover_sizes_dic
    """
    cover_sizes_key = misc.gen_key(cover_sizes)
    if cover_sizes_key not in ws_ex_cover_sizes_dic or len(ws_ex_cover_sizes_dic[cover_sizes_key]) < tuple_num:
        # entry not available or not enough samples. compute ws_ex_covers
        ws_num = len(ws_k_cover_dic[cover_sizes[0]])  # number of weighted samplings
        tsize = len(cover_sizes) # number of genes
        random_indices = choose_random_tuples(ws_num, len(cover_sizes), tuple_num)  # choose random indices
        ws_k_covers = [ws_k_cover_dic[k] for k in cover_sizes_key]  # extract weighted sampling for k's
        random_covers = \
            [[ws_k_covers[i][rtuple[i]] for i in range(tsize)] for rtuple in random_indices] # NULL covers
        ws_ex_cover_sizes = \
            [len(compute_ex_cover(random_covers[i])) for i in range(len(random_covers))] # NULL ex_cover_sizes
        ws_ex_cover_sizes.sort() # for efficiency
        ws_ex_cover_sizes_dic[cover_sizes_key] = ws_ex_cover_sizes # update ws_ex_cover_sizes_dic
    else:  # reuse if already in the dictionary
        ws_ex_cover_sizes = ws_ex_cover_sizes_dic[cover_sizes_key]
    return ws_ex_cover_sizes, ws_ex_cover_sizes_dic


# todo: check the parameters when call this function
def compute_me_co_pv_ws(covers, ws_k_cover_dic, me_co="me",
                     min_rank=100, init_pair_num=10**3, max_pair_num=10**6, ws_ex_cover_sizes_dic={}):
    """
    compute pvalue of ME/CO based on weighted sampling
    can compute pv for an arbitrary number of genes (not just a pair)
    :param covers [cover1, cover2, ..] : original over sizes list for genes to compute pv
    :param ws_k_cover_dic: k -> list of covers
    :param me_co: me or co
    :param min_rank: minimum rank (precision)
    :param init_pair_num: initial number of random pairs to start
    :param max_pair_num: max number of random pairs
    :param ws_ex_cover_sizes_dic: to reuse the previous sampling
                (k tuples) -> list of ex cover sizes
    :return ws_ps: pvalue
    :return ws_ex_cover_sizes_dic
    """

    cover_sizes = [len(c) for c in covers]
    # original cover size
    ex_cover_size = len(compute_ex_cover(covers))
    # compute random cover size (initial sampling)
    ws_cover_sizes, ws_ex_cover_sizes_dic = \
        compute_ws_cover_sizes(cover_sizes, ws_k_cover_dic, init_pair_num, ws_ex_cover_sizes_dic)
    while True:  # dynamic sampling
        # compute pvalue
        ws_rank = compute_rank(ex_cover_size, ws_cover_sizes, me_co, ordered=True)
        pair_num = len(ws_cover_sizes)
        if ws_rank < min_rank and pair_num < max_pair_num:  # if not precise enough
            # sample more pairs for better precision
            ws_cover_sizes, ws_ex_cover_sizes_dic = \
                compute_ws_cover_sizes(cover_sizes, ws_k_cover_dic, pair_num*10, ws_ex_cover_sizes_dic)
        else:
            ws_pv = ws_rank/float(pair_num)
            break
    return ws_pv, ws_ex_cover_sizes_dic


# todo: delete??
def compute_me_pvalue(ws_covers, cover, ordered=False, mid_pv=True):
    """
    given random instances (weighted sampling) and the origial cover size,
    compute empirical pvalue for mutual exclusivity
    compute mid pvalue if mid_pv is True
    mid pvalue is avg of strictly greater and geater/eqaul

    :param ws_covers: list of random cover sizes
    :param cover: original cover size
    :param ordered: ws_covers is sorted if True
    :param mid_pv: mid pvalue is computed if True
    :return: pvalue: pvalue
    """
    if ordered:  # if ws_covers is sorted
        less = bisect.bisect_left(ws_covers, cover)
        if mid_pv is True:
            lesseq = bisect.bisect_right(ws_covers, cover)
            pvalue = 1-(less+lesseq)/(2.0*len(ws_covers))
        else:
            pvalue = 1-less/float(len(ws_covers))
    else:  # if ws_covers is not sorted
        greatereq = len(filter(lambda c: c >= cover, ws_covers))
        if mid_pv is True:
            greater = len(filter(lambda c: c > cover, ws_covers))
            pvalue = (greater+greatereq)/(2.0*len(ws_covers))
        else:
            pvalue = greatereq/float(len(ws_covers))
    return pvalue


# todo: delete??
def compute_co_pvalue(ws_covers, cover, ordered=False, mid_pv=True):
    """
    given random instances (weighted sampling) and the origial cover size,
    compute empirical pvalue for co-occurrence
    compute mid pvalue if mid_pv is True
    mid pvalue is avg of strictly greater and geater/eqaul

    :param ws_covers: list of random cover sizes
    :param cover: original cover size
    :param ordered: ws_covers is sorted if True
    :param mid_pv: mid pvalue is computed if True
    :return: pvalue: pvalue
    """
    if ordered:  # if ws_covers is sorted
        lesseq = bisect.bisect_right(ws_covers, cover)
        if mid_pv is True:
            less = bisect.bisect_left(ws_covers, cover)
            pvalue = (less+lesseq)/(2.0*len(ws_covers))
        else:
            pvalue = lesseq/float(len(ws_covers))
    else:  # if ws_covers is not sorted
        lesseq = len(filter(lambda c: c <= cover, ws_covers))
        if mid_pv is True:
            less = len(filter(lambda c: c < cover, ws_covers))
            pvalue = (less+lesseq)/(2.0*len(ws_covers))
        else:
            pvalue = lesseq/float(len(ws_covers))
    return pvalue

