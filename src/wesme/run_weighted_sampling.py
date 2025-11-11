"""
generate random weighted sampling for each cancer/mutation type

>> python run_weighted_sampling.py cantype mtype rnum (-all)

  e.g.  to generate 1000 random samplings for BRCA somatic mutations
        only for observed cover sizes

>> python run_weighted_sampling.py BRCA smut 1000

  e.g.  to generate for all k [1... max_k]

>> python run_weighted_sampling.py BRCA smut 1000 -all

"""

import argparse
import logging
import numpy as np
from pathlib import Path
from data_proc import io_utils
from config import wrs_dir, data_dir  # directory to store weighted random sampling
logging.basicConfig(level=logging.INFO)


def comp_all_ks(mut_dic):
    """
    return all observed k's in the mutation profile
    k is the number of mutated samples for any gene
    :param mut_dic: gene -> number of mutated samples
    :return: set of k's
    """
    return set([len(x) for x in mut_dic.values()])

parser = argparse.ArgumentParser()
parser.add_argument("can", help="cancer type", type=str)
parser.add_argument("mtype", help="mutation type (smut, amp, del, cnv or alt)", type=str)
parser.add_argument("rnum", help="number of random instances to generate", type=int)
parser.add_argument("-m", "--mut_input", type=str, help="mutation file name")
parser.add_argument("-f", "--freq_input", type=str, help="frequency file name")
parser.add_argument("-o", "--output_dir", type=str, help="output directory")
parser.add_argument("--all", help="generate for all possible k's", action='store_true')
args = parser.parse_args()
can, mtype, rnum = args.can, args.mtype, args.rnum

# INPUT_FILES
if args.mut_input is not None:
    mut_file_name = args.mut_input
else:
    mut_file_name = data_dir+can+"/" + "_".join([can, mtype]) + "_list.txt"
if args.freq_input is not None:
    freq_file_name = args.freq_input
else:
    freq_file_name = data_dir+can+"/" + "_".join([can, mtype]) + "_freq.txt"

logging.info("%s", can)

# OUTPUT FILES/DIRECTORIES
# file name: wr_k_rnum.txt for each k
# file format: each line is a list comma separated sample indices
# obtained using weighted random sampling without replacement
if args.output_dir is not None:
    cur_wrs_dir = args.output_dir
else:
    # cur_wrs_dir = wrs_dir + can + "/" + mtype + "/"
    cur_wrs_dir = wrs_dir + can + "/"

# if output directory does not exist, create it
if not Path(cur_wrs_dir).exists():
    Path(cur_wrs_dir).mkdir(parents=True, exist_ok=True)

# read files
mut_list_dic, samples = io_utils.read_mut_list(mut_file_name, samples=True)
freq_dic = io_utils.read_dic(freq_file_name)
freqs = [float(freq_dic[s]) for s in samples]  # freq list

nsamples = len(samples)
if args.all: # for all i in [1..k]
    all_ks = range(nsamples)
    logging.info("for all possible [1-k]")
else: # only for observed k's
    all_ks = comp_all_ks(mut_list_dic)
    logging.info("for all observed k's")

# generate weighted random sampling for each k for rnum times
# and write to output files
for k in all_ks:
    if k == 0:  # no need to sample for k=0
        continue
    ws_file = open(cur_wrs_dir + "_".join(["wr", str(k), str(rnum)])+".txt", 'w')
    for i in range(rnum):
        wsamples = np.random.choice(range(nsamples), k, p=freqs, replace=False)
        ws_file.write("%s\n" % ",".join([str(idx) for idx in wsamples]))
    ws_file.close()

