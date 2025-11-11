"""
miscellaneous utility functions
"""


def get_positives(in_list):
    """ given a list of values, return positive indices
    """
    return list(filter(lambda i: in_list[i] > 0, range(len(in_list))))


def gen_key(alist):
    """
    sort the list and return in a tuple to use as a key
    :type alist: list
    :param alist: list of values
    :return: return sorted list as a tuple
    """
    alist = list(alist)
    alist.sort()
    return tuple(alist)


def gen_uniq_key(alist):
    """
    sort the list and return in a tuple to use as a key
    remove if there are redundancy (elements should be unique)
    :type alist: list
    :param alist: list of values
    :return: return sorted list as a tuple
    """
    alist = list(set(alist)) # remove redundancy
    alist.sort()
    return tuple(alist)
