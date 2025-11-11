"""
transform mutation matrix to mutation list

>> python mut_mat_to_list.py cantype mtype

    e.g.  to run for BRCA somatic mutations matrix
    
>> python mut_mat_to_list.py BRCA smut

"""
from pathlib import Path
import argparse
import logging
import pandas as pd
import data_proc.io_utils as io

logging.basicConfig(level=logging.DEBUG)

# read arguments
parser = argparse.ArgumentParser()
parser.add_argument("can", help="cancer type", type=str)
parser.add_argument("mtype", help="mutation type (smut, amp, del, cnv or alt)", type=str)
parser.add_argument("-i", "--input", type=str, help="input file name")  
parser.add_argument("-o", "--output", type=str, help="output file name")

args = parser.parse_args()
can, mtype = args.can, args.mtype

# Directory
data_dir = Path(__file__).resolve().parent.parent.parent / "data" / "processed"/ can
logging.debug("data directory: %s" %data_dir)

# INPUT_FILES
if args.input is not None:
    mat_file = args.input
else:
    mat_file = data_dir / "mut_data.tsv"
    
# OUTPUT_FILES
# each line: gene\tsample1, sample2, ... (first line: sample\tsample_name1, sample_name2, ...)
if args.output is not None:
    mut_file = args.output
else:
    mut_file = "data/" + can + "/" + "_".join([can, mtype]) + "_list.txt"
    
# if output directory does not exist, create it
if not Path(mut_file).parent.exists():
    Path(mut_file).parent.mkdir(parents=True, exist_ok=True)
    
logging.debug("mutation file: %s" %mut_file)

# read mutation matrix
mut_mat = pd.read_csv(mat_file, sep="\t", index_col=0)

# extract samples and mutation dictionary
samples = mut_mat.columns.tolist()
mut_dic = {gene: mut_mat.loc[gene].tolist() for gene in mut_mat.index}

# write mutation list
io.write_mut_list(mut_dic, mut_file, samples=samples)