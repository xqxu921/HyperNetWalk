"""
compute FDR using permutation and binning

>> python cantype mtype me_co pnum list_of_bin_ths

e.g., compute FDR for each ME pvalues using 100 permuted instances
    with 2 gene groups (= 3 bins for pairs)
    mr_th = 2% of samples

>> python comp_fdr.py BRCA smut me 100 2

"""
import os
from pathlib import Path
import pandas
import logging
import argparse
import math
logging.basicConfig(level=logging.INFO)


from config import pre_dir, results_dir
from data_proc import io_utils

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("can", help="cancer type", type=str)
parser.add_argument("mtype", help="mutation type", type=str)
parser.add_argument("me_co", help="me or co", type=str)
parser.add_argument("pnum", help="number of permuted instances to be used for FDR computation", type=int)
parser.add_argument("bin_ths", help="percentage mutation rates to divide bins", nargs='*')
parser.add_argument("-m", "--mut_input", type=str, help="mutation file name")
parser.add_argument("-p", "--pv_input", type=str, help="pv file name")
parser.add_argument("-n", "--null_prefix", type=str, help="null file name prefix")
parser.add_argument("-f", "--fdr_output", type=str, help="fdr output file name")
parser.add_argument("-pf", "--pv_fdr_output", type=str, help="pv/fdr output file name")

args = parser.parse_args()
cantype, mtype, me_co, pnum = args.can, args.mtype, args.me_co, args.pnum
bin_ths = [float(x) for x in args.bin_ths]
logging.info("comp %s for %s, %s (bin_th=%s)" %(me_co, cantype, mtype, ",".join(args.bin_ths)))

# input files
if args.mut_input is not None:
    mut_file = args.mut_input
else:
    mut_file = "data/" + cantype + "/" + cantype + "_" + mtype + "_list.txt"

if args.pv_input is not None:
    orig_file = args.pv_input
else:
    orig_file = "results/" + cantype + "/" + "_".join([cantype, mtype, me_co, "pvs_0.1.txt"])

if args.null_prefix is not None:
    null_pv_file_prefix = args.null_prefix
else:
    null_results_dir = pre_dir + "permuted_pv/"+cantype+"/"
    null_pv_file_prefix = null_results_dir + "_".join([cantype, mtype,"perm", "pv_"])

# output files
prefix_list = [cantype, mtype, me_co]
suffix_list = [str(len(bin_ths)+1)+"bins"]+args.bin_ths # number of bins followed by each threshold
if args.fdr_output is not None:
    fdr_file = args.fdr_output
else:
    fdr_file = results_dir + cantype + "/" + "fdr/" + "_".join(prefix_list+["fdr"]+suffix_list) + ".txt"
if args.pv_fdr_output is not None:
    pv_fdr_file = args.pv_fdr_output
else:
    pv_fdr_file = results_dir + cantype + "/" + "fdr/" + "_".join(prefix_list+["pv_fdr"]+suffix_list) + ".txt"

# if the output directory does not exist, create it
if not os.path.exists(os.path.dirname(fdr_file)):
    os.makedirs(os.path.dirname(fdr_file))
if not os.path.exists(os.path.dirname(pv_fdr_file)):
    os.makedirs(os.path.dirname(pv_fdr_file))

def add_mr(df, mr_dic):
    """
    add mutation rate information an additional columns in the original dataframe
    :param df (data frame): pvalues
    :param mr_dir: muration rate dic
    :return: modified df
    """
    df["mr1"] = df["gene1"].map(lambda g: mr_dic[g])
    df["mr2"] = df["gene2"].map(lambda g: mr_dic[g])
    return df


def select_pvs_mr(df, mr1, mr2):
    """
    select pairs with mutation rate falls in the condition given by mr1, mr2
    :param df (data frame): pvalues
    :param mr1: (start, end) mutation rate for gene 1
    :param mr2: (start, end) mutation rate for gene 2
    :return: selected pvalues,  df for selected genes
    """
    selected_df = df[((df['mr1'] >= mr1[0]) & (df['mr2'] >=mr2[0]) & ((df['mr1'] < mr1[1]) & (df['mr2'] < mr2[1])))|
                          ((df['mr1'] >= mr2[0]) & (df['mr2'] >=mr1[0]) & ((df['mr1'] < mr2[1]) & (df['mr2'] < mr1[1])))]
    selected_pvs = selected_df["pv"].tolist()
    selected_pvs.sort()
    return selected_pvs, selected_df


def comp_fdr(orig_pvs, rand_pvs, pnum):
    """
    compute fdr
    :param orig_pvs: original pvalues for selected pairs
    :param rand_pvs: NULL pvlaues for the pairs in the same bin
    :return: fdr: pv->fdr dic
    """
    fdr = {}
    # add dummy
    rand_pvs.append(10)
    j = 0
    pnum = float(pnum)
    for i in range(len(orig_pvs)):
        pv = orig_pvs[i]
        # count number of random pvs that are less than or equal to pv
        # move j from left to right
        while rand_pvs[j] <= pv:
            j += 1
        fdr[pv] = j/((i+1)*pnum)
    return fdr


# compute mutation rates and find highly mutated genes
mut_list_dic, samples = io_utils.read_mut_list(mut_file, samples=True)
nsamples = len(samples)
logging.info("%d samples" %nsamples)

# mutation rates information
mr_genes_dic = {}
mr_dic = dict([(gene, len(mut_list_dic[gene])) for gene in mut_list_dic])

# read original pvalues
orig = pandas.read_table(orig_file)
orig_df = add_mr(orig, mr_dic)
orig_df = orig_df.rename(columns={"pv (ws)":"pv"})

# read NULL pvalues
all_rands = []
for i in range(pnum):
    rand_file = null_pv_file_prefix + str(i) + ".txt"
    rand = pandas.read_table(rand_file)
    rand = rand.rename(columns={"pv (ws)":"pv"})
    all_rands.append(add_mr(rand, mr_dic))

mr_ths = [0.01*x for x in bin_ths]
mrs = [0]+mr_ths+[1] # add extra dummy for binning range
mnum = [math.ceil(nsamples*mr) for mr in mrs]
logging.info("%s\n" %str(mrs))

# to fix the memory problem
# rand_df = pandas.concat(all_rands)
fdr_dfs = []
orig_df_fdrs = []
selected_rand_df = [[] for k in range(pnum)]
bin_info = []
# compute fdr using binning
for i in range(len(mrs)-1):
    mr1 = (mnum[i], mnum[i+1])  # range of smaller mr
    for j in range(i, len(mrs)-1):
        mr2 = (mnum[j], mnum[j+1])  # range of bigger mr
        orig_pvs, selected_orig_df = select_pvs_mr(orig_df, mr1, mr2)  # selected gene pair in mr1, mr2
        if len(orig_pvs) == 0: # if the bin is empty
            continue
        # select NULL pvs in mr1, mr2
        rand_pvs_all = []
        for k in range(pnum):
            temp_pvs, selected_df = select_pvs_mr(all_rands[k], mr1, mr2)
            rand_pvs_all.extend(temp_pvs)
        rand_pvs_all.sort()
        fdr = comp_fdr(orig_pvs, rand_pvs_all, pnum)
        pvs = fdr.keys()
        pvs = list(pvs)
        pvs.sort()
        # make a data frame with pv and fdr (for column concatenation later)
        rows = [(pv, fdr[pv]) for pv in pvs]
        mr_label = "mr="+str((mrs[i], mrs[i+1]))+","+str((mrs[j], mrs[j+1]))
        fdr_df = pandas.DataFrame(data=rows, columns=["pv", mr_label])
        fdr_dfs.append(fdr_df)
        # make a data frame with pv and fdr (for row concatenation)
        orig_df_fdrs.append(pandas.merge(selected_orig_df, fdr_df, how='left', on=['pv']).rename(columns={mr_label:"fdr"}))
        bin_info.append((mr1,mr2))

# create two data frames and 2 output files
# merge by columns (pv -> fdr mapping)
fdr_df = fdr_dfs[0]
for df in fdr_dfs[1:]:
    fdr_df = fdr_df.merge(df, how='outer')

# merge by rows (original pv dataframe + fdr)
orig_df_fdr = pandas.concat(orig_df_fdrs)

# write the results
sorted_fdr = fdr_df.sort_values(["pv"])
sorted_fdr.to_csv(fdr_file, sep="\t", index=False)

sorted_df_fdr = orig_df_fdr.sort_values(["fdr"])
fdr_th=0.5 # filter with FDR_th and write the results
sorted_df_fdr[sorted_df_fdr['fdr'] <= fdr_th].to_csv(pv_fdr_file, sep="\t", index=False)

# bin_info 和 orig_df_fdrs 顺序一致
# bin_info = [(mr1, mr2), ...]
# orig_df_fdrs = [df1, df2, ...]

# 保存 bin_info 到文件
with open(results_dir + cantype + "/" + "fdr/" + "_".join(prefix_list+["bin_info"]+suffix_list) + ".txt", "w") as f:
    f.write("mr1_start\tmr1_end\tmr2_start\tmr2_end\n")
    for mr1, mr2 in bin_info:
        f.write(f"{mr1[0]}\t{mr1[1]}\t{mr2[0]}\t{mr2[1]}\n")
# 找到 (R,H) 和 (H,H) 的索引
rh_idx = None
hh_idx = None
for idx, (mr1, mr2) in enumerate(bin_info):
    if mr1 == (mnum[0], mnum[1]) and mr2 == (mnum[1], mnum[2]):
        rh_idx = idx
    elif mr1 == (mnum[1], mnum[2]) and mr2 == (mnum[1], mnum[2]):
        hh_idx = idx

# 取出对应的 DataFrame
rh_df_fdr = orig_df_fdrs[rh_idx] if rh_idx is not None else pandas.DataFrame()
hh_df_fdr = orig_df_fdrs[hh_idx] if hh_idx is not None else pandas.DataFrame()

# 后续筛选和保存
filtered_rh_df_fdr = rh_df_fdr[(rh_df_fdr['pv']<=0.001) & (rh_df_fdr['fdr']<=0.125)] if not rh_df_fdr.empty else pandas.DataFrame()
filtered_hh_df_fdr = hh_df_fdr[(hh_df_fdr['pv']<=0.01) & (hh_df_fdr['fdr']<=0.125)] if not hh_df_fdr.empty else pandas.DataFrame()

filtered_orig_df_fdr = pandas.concat([filtered_rh_df_fdr, filtered_hh_df_fdr])
filtered_orig_df_fdr.to_csv(results_dir + cantype + "/" + "fdr/" + "_".join(prefix_list+["filtered_pv_fdr"]+suffix_list) + ".txt", sep="\t", index=False)
filtered_orig_df_fdr.to_csv(Path(__file__).resolve().parent.parent.parent / "data" / "NETWORK" / "_".join(cantype+["me"]+["net"])+".txt",sep="\t", index=False)
logging.info("done")
        
# # (R,H) 
# rh_df_fdr = orig_df_fdrs[0]
# # 保留pv<0.001且fdr<0.125的行
# filtered_rh_df_fdr = rh_df_fdr[(rh_df_fdr['pv']<=0.001) & (rh_df_fdr['fdr']<=0.125)]

# # (H,H)
# hh_df_fdr = orig_df_fdrs[1]
# # 保留pv<0.01的行
# filtered_hh_df_fdr = hh_df_fdr[(hh_df_fdr['pv']<=0.01) & (hh_df_fdr['fdr']<=0.125)]

# filtered_orig_df_fdr = pandas.concat([filtered_rh_df_fdr, filtered_hh_df_fdr])
# filtered_orig_df_fdr.to_csv(results_dir + cantype + "/" + "fdr/" + "_".join(prefix_list+["filtered_pv_fdr"]+suffix_list) + ".txt", sep="\t", index=False)
# logging.info("done")