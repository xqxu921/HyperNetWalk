"""
compute mutual exclusivity rank/pv via misc methods
"""

import scipy.stats as ss


def compute_jaccard(cover1, cover2):
    """
    compute jaccard index given two sets
    :param cover1:
    :param cover2:
    :return: ji
    """
    union_cov = set(cover1).union(cover2)
    cross_cov = set(cover1).intersection(cover2)
    return len(cross_cov)/float(len(union_cov))


def compute_me_pv_hypergeom(cover1, cover2, nsamples):
    """
    compute ME pvalue using hypergeometric test
    :param cover1:
    :param cover2:
    :param nsamples: total number of samples
    :return: h_pv
    """
    param = (len(set(cover1).intersection(cover2)), nsamples, len(cover1), len(cover2))
    h_pv = ss.hypergeom.cdf(*param)
    return h_pv


def compute_co_pv_hypergeom(cover1, cover2, nsamples):
    """
    compute CO pvalue using hypergeometric test
    :param cover1:
    :param cover2:
    :param nsamples: total number of samples
    :return: h_pv
    """
    param = (len(cover2)-len(set(cover1).intersection(cover2)), nsamples, nsamples-len(cover1), len(cover2))
    h_pv = ss.hypergeom.cdf(*param)
    return h_pv


def compute_hypergeom(cover1, cover2, nsamples):
    """
    compute "PMF" with hypergometric
    :param cover1:
    :param cover2:
    :param nsamples: total number of samples
    :return: h_pv
    """
    param = (len(set(cover1).intersection(cover2)), nsamples, len(cover1), len(cover2))
    h_pv = ss.hypergeom.pmf(*param)
    return h_pv
