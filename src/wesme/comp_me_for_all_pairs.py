"""
compute WeSME pvalue for all pairs
each mutation type separately

>> python comp_me_for_all_pairs.py BRCA smut 0.1

"""

import argparse
import os
import logging
import pandas
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

from config import wrs_dir
from mut_ex import ws_mut_ex, misc_mut_ex
from data_proc import io_utils


# arguments
parser = argparse.ArgumentParser()
parser.add_argument("can", help="cancer type", type=str)
parser.add_argument("mtype", help="mutation type (smut, amp, del, cnv or alt)", type=str)
parser.add_argument("pth", help="pvalue threshold to write in the output file", type=float)
parser.add_argument("-m", "--mut_input", type=str, help="mutation file name")
parser.add_argument("-s", "--sampling_dir", type=str, help="sampling directory")
parser.add_argument("-o", "--output", type=str, help="output file name")

args = parser.parse_args()
cantype, mtype, pth = args.can, args.mtype, args.pth
logging.debug("comp ME for %s, %s" %(cantype, mtype))

# INPUT_FILES/DIRECTORIES
if args.mut_input is not None:
    mut_file = args.mut_input
else:
    mut_dir = "data/" + cantype + "/"
    mut_file = mut_dir + "_".join([cantype, mtype]) + "_list.txt"
if args.sampling_dir is not None:
    cur_wrs_dir = args.sampling_dir
else:
    cur_wrs_dir = wrs_dir + cantype + "/"

# OUTPUT_FILES
if args.output is not None:
    pv_file = args.output
else:
    results_dir = "results/" + cantype + "/"
    pv_file = results_dir + cantype + "_" + mtype + "_me_pvs_"+str(pth)+".txt"

# if output directory does not exist, create it
if not os.path.exists(os.path.dirname(pv_file)):
    os.makedirs(os.path.dirname(pv_file))

# read mutation profile
alt_list_dic, samples = io_utils.read_mut_list(mut_file, samples=True)
gene_ks = dict([(gene, len(alt_list_dic[gene])) for gene in alt_list_dic])  # gene -> cover size dic
nsamples = len(samples)

# read weighted sampling files for all k
# use the weighted sampling files in cur_wrs_dir
kfiles = os.listdir(cur_wrs_dir)
# store weighted sampling for each k: k --> list of arrays
ws_k_cover_dic = {}
for kf in kfiles:
    # read random file for k
    k = int(kf.split(".")[0].split("_")[1])
    ws_k_cover_dic[k], samples = io_utils.read_mut_list(cur_wrs_dir + kf, genes=False)

ws_pv_dic = {}  # WeSME
jaccard_dic = {}  # jaccard index

# compute ME for all pairs
genes = list(alt_list_dic.keys())
ws_ex_cover_sizes_dic = {}
for i in range(len(genes)):
    gene1 = genes[i]
    if i%1000 == 0:
        logging.info("%d " % i)
    for j in range(i+1, len(genes)):
        gene2 = genes[j]
        gene_pairs = (gene1, gene2)
        # compute pvalue using weighted sampling
        covers = [alt_list_dic[gene1], alt_list_dic[gene2]]
        ws_pv, ws_ex_cover_sizes_dic  = ws_mut_ex.compute_me_co_pv_ws\
            (covers, ws_k_cover_dic, max_pair_num=10**5, ws_ex_cover_sizes_dic=ws_ex_cover_sizes_dic)
        if ws_pv > pth:  # if pvalue is not significant enough
            continue
        ws_pv_dic[gene_pairs] = ws_pv # WeSME pvalue
        # compute jaccard index
        jaccard_dic[gene_pairs] = misc_mut_ex.compute_jaccard(alt_list_dic[gene1], alt_list_dic[gene2])


# write the results
rows = []
for gp in ws_pv_dic:
    rows.append((gp[0], gp[1], jaccard_dic[gp], ws_pv_dic[gp]))

# write the results
all_pvs = pandas.DataFrame(data=rows, columns=["gene1", "gene2", "jaccard index", "pv (ws)"])
all_pvs.to_csv(pv_file, sep="\t", index=False)
# """
# compute WeSME pvalue for all pairs
# each mutation type separately

# >> python comp_me_for_all_pairs.py BRCA smut 0.1 -n 4

# 智能并行策略：
# - 按(k1,k2)组合分组基因对
# - 每个组在同一进程内处理，最大化缓存复用
# - 理论上能接近顺序执行的缓存命中率，同时获得并行加速
# """

# import argparse
# import os
# import logging
# import pandas
# import time
# import multiprocessing as mp
# from tqdm import tqdm
# from collections import defaultdict

# logging.basicConfig(
#     level=logging.INFO,
#     format='%(asctime)s - %(levelname)s - %(message)s'
# )
# logger = logging.getLogger()

# from config import wrs_dir
# from mut_ex import ws_mut_ex, misc_mut_ex
# from data_proc import io_utils


# def compute_gene_pairs_for_k_group(args):
#     """
#     计算同一(k1,k2)组合的所有基因对
#     这个函数内的所有计算会共享缓存！
#     """
#     k_pair, gene_pairs, alt_list_dic, ws_k_cover_dic, pth = args
    
#     results = []
#     ws_ex_cover_sizes_dic = {}  # 这个组内共享的缓存
    
#     cache_hits = 0
#     cache_misses = 0
    
#     for gene1, gene2 in gene_pairs:
#         covers = [alt_list_dic[gene1], alt_list_dic[gene2]]
#         cover_sizes_key = tuple(sorted([len(c) for c in covers]))
        
#         # 追踪缓存
#         if cover_sizes_key in ws_ex_cover_sizes_dic:
#             cache_hits += 1
#         else:
#             cache_misses += 1
        
#         ws_pv, ws_ex_cover_sizes_dic = ws_mut_ex.compute_me_co_pv_ws(
#             covers, ws_k_cover_dic, max_pair_num=10**5,
#             ws_ex_cover_sizes_dic=ws_ex_cover_sizes_dic
#         )
        
#         if ws_pv <= pth:
#             jaccard = misc_mut_ex.compute_jaccard(
#                 alt_list_dic[gene1], alt_list_dic[gene2]
#             )
#             results.append((gene1, gene2, jaccard, ws_pv))
    
#     # 返回统计信息
#     total = cache_hits + cache_misses
#     hit_rate = cache_hits / total * 100 if total > 0 else 0
    
#     return {
#         'k_pair': k_pair,
#         'results': results,
#         'total_pairs': len(gene_pairs),
#         'cache_hits': cache_hits,
#         'cache_misses': cache_misses,
#         'hit_rate': hit_rate,
#         'unique_keys': len(ws_ex_cover_sizes_dic)
#     }


# def group_gene_pairs_by_coverage(alt_list_dic):
#     """
#     按(k1,k2)覆盖度组合对基因对进行分组
#     返回: {(k1, k2): [(gene1, gene2), ...]}
#     """
#     # 首先按覆盖度对基因分类
#     coverage_to_genes = defaultdict(list)
#     for gene, samples in alt_list_dic.items():
#         k = len(samples)
#         coverage_to_genes[k].append(gene)
    
#     # 生成所有(k1,k2)组合的基因对
#     k_values = sorted(coverage_to_genes.keys())
#     k_pair_to_gene_pairs = defaultdict(list)
    
#     for i, k1 in enumerate(k_values):
#         genes1 = coverage_to_genes[k1]
        
#         # 同一覆盖度的基因对
#         for idx1 in range(len(genes1)):
#             for idx2 in range(idx1 + 1, len(genes1)):
#                 k_pair = (k1, k1)
#                 k_pair_to_gene_pairs[k_pair].append((genes1[idx1], genes1[idx2]))
        
#         # 不同覆盖度的基因对
#         for k2 in k_values[i + 1:]:
#             genes2 = coverage_to_genes[k2]
#             k_pair = (k1, k2)
#             for g1 in genes1:
#                 for g2 in genes2:
#                     k_pair_to_gene_pairs[k_pair].append((g1, g2))
    
#     return k_pair_to_gene_pairs


# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument("can", help="cancer type", type=str)
#     parser.add_argument("mtype", help="mutation type (smut, amp, del, cnv or alt)", type=str)
#     parser.add_argument("pth", help="pvalue threshold to write in the output file", type=float)
#     parser.add_argument("-m", "--mut_input", type=str, help="mutation file name")
#     parser.add_argument("-s", "--sampling_dir", type=str, help="sampling directory")
#     parser.add_argument("-o", "--output", type=str, help="output file name")
#     parser.add_argument("-n", "--n_cores", type=int, default=1,
#                        help="number of cores (1=sequential, >1=parallel)")
#     parser.add_argument("-v", "--verbose", action="store_true", help="show detailed statistics")

#     args = parser.parse_args()
#     cantype, mtype, pth = args.can, args.mtype, args.pth
#     n_cores = args.n_cores
    
#     start_time = time.time()
#     logger.info(f"Computing ME for {cantype}, {mtype}, p-value threshold={pth}")

#     # INPUT
#     if args.mut_input is not None:
#         mut_file = args.mut_input
#     else:
#         mut_dir = "data/" + cantype + "/"
#         mut_file = mut_dir + "_".join([cantype, mtype]) + "_list.txt"
    
#     if args.sampling_dir is not None:
#         cur_wrs_dir = args.sampling_dir
#     else:
#         cur_wrs_dir = wrs_dir + cantype + "/"

#     # OUTPUT
#     if args.output is not None:
#         pv_file = args.output
#     else:
#         results_dir = "results/" + cantype + "/"
#         pv_file = results_dir + cantype + "_" + mtype + "_me_pvs_" + str(pth) + ".txt"

#     if not os.path.exists(os.path.dirname(pv_file)):
#         os.makedirs(os.path.dirname(pv_file))

#     # Read data
#     logger.info("Reading mutation data...")
#     alt_list_dic, samples = io_utils.read_mut_list(mut_file, samples=True)
#     nsamples = len(samples)
#     logger.info(f"Loaded {len(alt_list_dic)} genes across {nsamples} samples")

#     logger.info("Reading weighted sampling data...")
#     kfiles = os.listdir(cur_wrs_dir)
#     ws_k_cover_dic = {}
#     for kf in kfiles:
#         k = int(kf.split(".")[0].split("_")[1])
#         ws_k_cover_dic[k], _ = io_utils.read_mut_list(cur_wrs_dir + kf, genes=False)
    
#     # 统计覆盖度分布
#     coverages = [len(samples) for samples in alt_list_dic.values()]
#     unique_k = len(set(coverages))
#     logger.info(f"Coverage range: [{min(coverages)}, {max(coverages)}]")
#     logger.info(f"Unique coverage values: {unique_k}")
#     logger.info(f"Average coverage: {sum(coverages)/len(coverages):.1f}")
    
#     # 按(k1,k2)分组
#     logger.info("Grouping gene pairs by coverage combinations...")
#     k_pair_to_gene_pairs = group_gene_pairs_by_coverage(alt_list_dic)
    
#     total_pairs = sum(len(pairs) for pairs in k_pair_to_gene_pairs.values())
#     num_k_groups = len(k_pair_to_gene_pairs)
    
#     logger.info(f"Total gene pairs: {total_pairs:,}")
#     logger.info(f"Unique (k1,k2) groups: {num_k_groups}")
#     logger.info(f"Average pairs per group: {total_pairs/num_k_groups:.1f}")
    
#     # 显示最大的几个组
#     if args.verbose:
#         sorted_groups = sorted(k_pair_to_gene_pairs.items(), 
#                               key=lambda x: len(x[1]), reverse=True)
#         logger.info("Top 5 largest (k1,k2) groups:")
#         for k_pair, pairs in sorted_groups[:5]:
#             logger.info(f"  {k_pair}: {len(pairs):,} pairs")
    
#     all_results = []
    
#     if n_cores <= 1:
#         # 顺序执行
#         logger.info("Running SEQUENTIAL version...")
        
#         with tqdm(total=num_k_groups, desc="Processing (k1,k2) groups", unit="group") as pbar:
#             for k_pair, gene_pairs in k_pair_to_gene_pairs.items():
#                 result = compute_gene_pairs_for_k_group(
#                     (k_pair, gene_pairs, alt_list_dic, ws_k_cover_dic, pth)
#                 )
#                 all_results.append(result)
                
#                 if args.verbose:
#                     pbar.set_postfix({
#                         'k_pair': f'{k_pair}',
#                         'pairs': len(gene_pairs),
#                         'sig': len(result['results']),
#                         'cache_hit': f"{result['hit_rate']:.1f}%"
#                     })
#                 pbar.update(1)
    
#     else:
#         # 并行执行 - 按(k1,k2)组分配任务
#         logger.info(f"Running PARALLEL version with {n_cores} cores...")
#         logger.info(f"Strategy: Each (k1,k2) group processed by single worker")
#         logger.info(f"Expected: High cache hit rate within each group")
        
#         # 准备任务
#         tasks = [
#             (k_pair, gene_pairs, alt_list_dic, ws_k_cover_dic, pth)
#             for k_pair, gene_pairs in k_pair_to_gene_pairs.items()
#         ]
        
#         # 按任务大小排序，大任务优先（负载均衡）
#         tasks.sort(key=lambda x: len(x[1]), reverse=True)
        
#         # 并行处理
#         with mp.Pool(processes=n_cores) as pool:
#             results_iter = pool.imap(compute_gene_pairs_for_k_group, tasks)
            
#             with tqdm(total=num_k_groups, desc="Processing (k1,k2) groups", unit="group") as pbar:
#                 for result in results_iter:
#                     all_results.append(result)
                    
#                     if args.verbose:
#                         pbar.set_postfix({
#                             'k_pair': f'{result["k_pair"]}',
#                             'sig': len(result['results']),
#                             'hit': f"{result['hit_rate']:.1f}%"
#                         })
#                     pbar.update(1)
    
#     # 合并结果
#     final_results = []
#     total_cache_hits = 0
#     total_cache_misses = 0
#     total_unique_keys = 0
    
#     for result in all_results:
#         final_results.extend(result['results'])
#         total_cache_hits += result['cache_hits']
#         total_cache_misses += result['cache_misses']
#         total_unique_keys += result['unique_keys']
    
#     # 统计
#     elapsed = time.time() - start_time
#     overall_hit_rate = total_cache_hits / (total_cache_hits + total_cache_misses) * 100
    
#     logger.info("=" * 70)
#     logger.info("Performance Statistics:")
#     logger.info(f"  Total pairs computed: {total_pairs:,}")
#     logger.info(f"  Significant pairs found: {len(final_results):,} ({100*len(final_results)/total_pairs:.2f}%)")
#     logger.info(f"  (k1,k2) groups processed: {num_k_groups}")
#     logger.info(f"  Overall cache hit rate: {overall_hit_rate:.2f}%")
#     logger.info(f"  Total cache hits: {total_cache_hits:,}")
#     logger.info(f"  Total cache misses: {total_cache_misses:,}")
#     logger.info(f"  Total unique keys computed: {total_unique_keys}")
#     logger.info(f"  Total time: {elapsed:.2f}s ({elapsed/60:.2f}min)")
#     logger.info(f"  Processing speed: {total_pairs/elapsed:.0f} pairs/sec")
    
#     if n_cores > 1:
#         sequential_estimate = elapsed * n_cores * 0.9  # 估计的顺序执行时间
#         speedup = sequential_estimate / elapsed
#         logger.info(f"  Estimated speedup: {speedup:.2f}x")
    
#     logger.info("=" * 70)
    
#     # 详细统计（如果启用verbose）
#     if args.verbose and all_results:
#         logger.info("\nPer-group statistics:")
#         sorted_results = sorted(all_results, key=lambda x: x['total_pairs'], reverse=True)
#         for i, result in enumerate(sorted_results[:10]):
#             logger.info(f"  Group {i+1}: {result['k_pair']} - "
#                        f"{result['total_pairs']:,} pairs, "
#                        f"{len(result['results'])} significant, "
#                        f"cache hit: {result['hit_rate']:.1f}%, "
#                        f"unique keys: {result['unique_keys']}")

#     # Write results
#     if final_results:
#         all_pvs = pandas.DataFrame(
#             data=final_results,
#             columns=["gene1", "gene2", "jaccard index", "pv (ws)"]
#         )
#         all_pvs = all_pvs.sort_values('pv (ws)')
#         all_pvs.to_csv(pv_file, sep="\t", index=False)
        
#         logger.info(f"\nResults saved to: {pv_file}")
#         logger.info(f"P-value range: [{all_pvs['pv (ws)'].min():.2e}, {all_pvs['pv (ws)'].max():.2e}]")
#         logger.info(f"Jaccard index range: [{all_pvs['jaccard index'].min():.3f}, "
#                    f"{all_pvs['jaccard index'].max():.3f}]")
#     else:
#         logger.info("No significant gene pairs found")
#         all_pvs = pandas.DataFrame(columns=["gene1", "gene2", "jaccard index", "pv (ws)"])
#         all_pvs.to_csv(pv_file, sep="\t", index=False)


# if __name__ == "__main__":
#     main()