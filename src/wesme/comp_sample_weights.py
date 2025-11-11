"""
compute the frequency of each sample
>> python comp_sample_weights.py cantype mtype

  e.g.  to run for BRCA somatic mutations

>> python comp_sample_weights.py BRCA smut

"""

import argparse
import logging
from data_proc import io_utils
logging.basicConfig(level=logging.INFO) 


def count_muts(mut_list_dic, samples):
    """
    count the number of mutated genes for each sample
    :param must_list_dic: gene -> list of mutated sample indices
    :param samples: list of samples
    :return: count_dic: sample -> # altered genes
    """
    nsamples = len(samples)
    counts = [0 for i in range(nsamples)]
    # compute mut count for each sample
    for gene in mut_list_dic:
        for i in mut_list_dic[gene]:
            counts[i] += 1
    # convert to freq
    counts_dic = dict(zip(*[samples, counts]))
    return counts_dic


def compute_sample_freq(count_dic):
    """
    compute relative frequency for each sample
    :param count_dic:sample -> # altered genes
    :return: freq_dic:sample -> # altered genes/# total
    """
    total = float(sum(count_dic.values()))
    freq_dic = dict([(sample, count_dic[sample]/total) for sample in count_dic])
    return freq_dic

# read argument
parser = argparse.ArgumentParser()
parser.add_argument("can", help="cancer type", type=str)
parser.add_argument("mtype", help="mutation type (smut, amp, or del)", type=str)
parser.add_argument("-i", "--input", type=str, help="input file name")
parser.add_argument("-o", "--output", type=str, help="output file name")

args = parser.parse_args()
can, mtype = args.can, args.mtype

# Directory
data_dir = "data/"+can+"/"

# INPUT_FILES
if args.input is not None:
    smut_file = args.input
else:
    smut_file = data_dir + "_".join([can, mtype, "list.txt"])

# OUTPUT_FILES
# each line: sample_label\tfreq
if args.output is not None:
    smut_freq_file = args.output
else:
    smut_freq_file = data_dir + "_".join([can, mtype]) + "_freq.txt"

logging.info("%s", can)

# read mut file
smut_list_dic, samples = io_utils.read_mut_list(smut_file, samples=True)

smut_count = count_muts(smut_list_dic, samples)
smut_freq = compute_sample_freq(smut_count)
io_utils.write_dic(smut_freq, smut_freq_file)
