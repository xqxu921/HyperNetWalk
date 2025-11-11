
"""
permute a given mutation profile

>> python run_permute_cover_each_type.py cantype fid pnum

  e.g.  to run 100 * |E| swaps to create a permuted mutational profile
  for BRCA somatic mutations and create "BRCA_smut_permuted_cover_1.txt"

>> python permute_cover_each_type.py BRCA 1 100

modify INPUT/OUTPUT FILES as well as config.py file to change options.

"""

import argparse
import logging
logging.basicConfig(level=logging.DEBUG)
import time
import os
from config import data_dir, permuted_dir
from mut_ex import permute_mut_data

import data_proc.io_utils as io

# read arguments
parser = argparse.ArgumentParser()
parser.add_argument("can", help="cancer type", type=str)
parser.add_argument("mtype", help="mutation type (smut, amp, del, cnv or alt)", type=str)
parser.add_argument("fid", help="file id", type=int)
parser.add_argument("pnum", help="number of permutations", type=int)
parser.add_argument("-i", "--input", type=str, help="input file name")
parser.add_argument("-p", "--output_prefix", type=str, help="output file name prefix")

args = parser.parse_args()
cantype, mtype, fid, pnum = args.can, args.mtype, args.fid, args.pnum,

logging.info("%d-th permutation.. (permute %d times)" % (fid, pnum))


# INPUT_FILES/DIRECTORIES

# INPUT_FILES
if args.input is not None:
    mut_file = args.input
else:
    mut_dir = data_dir + cantype + "/"
    mut_file = mut_dir + "_".join([cantype, mtype]) + "_list.txt"

# OUTPUT_FILES
if args.output_prefix is not None:
    output_file = args.output_prefix + "_"+str(fid)+".txt"
else:
    # cur_permuted_dir = permuted_dir + cantype + "/" + mtype + "/"
    cur_permuted_dir = permuted_dir + cantype + "/"
    output_file = cur_permuted_dir+"_".join([cantype, mtype])+"_permuted_cover_"+str(fid)+".txt"

# if output directory does not exist, create it
if not os.path.exists(os.path.dirname(output_file)):
    os.makedirs(os.path.dirname(output_file))

# read data
logging.info("\n****** reading mutation data files (%s)...\n" %cantype)
mut_list_dic, samples = io.read_mut_list(mut_file, samples=True)
nsamples = len(samples)
genes = mut_list_dic.keys()

# todo: simplify to be run for one cancer type?
# functions designed to be run for multiple cancer types
# run it only for one cancer type
type_idx_dic = dict([(cantype, range(nsamples))])
cancers = [cantype]

# construct a bipartite graph for each cancer type
mut_graphs = permute_mut_data.construct_mut_graph_per_type_from_list(mut_list_dic, cancers, type_idx_dic)

# permute each bipartite graph separately
permuted_graphs = {}
permuted_graphs[cantype] = permute_mut_data.permute_mut_graph(mut_graphs[cantype], genes, range(nsamples), pnum)

# merge muted graphs and construct the permuted dic
permuted_mut_dic = permute_mut_data.construct_mut_dic_from_graphs(permuted_graphs, genes, nsamples)

io.write_mut_list(permuted_mut_dic, output_file, samples=samples)
