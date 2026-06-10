"""
说明：main_modified.py 相对于 main_adjusted.py 的主要差异如下：

1. ✅ mutation matrix 的基因过滤标准不同：
   - main_adjusted.py：仅保留突变频率 ≥ 2% 的基因（即至少出现在 2% 样本中的突变基因）；
     ```python
     mutationFrequency = (mutMatrix.values != 0).sum(axis=1) / numSamples
     mutMatrix = mutMatrix.loc[mutationFrequency >= 0.02]
     ```
   - main_modified.py：保留所有在任意样本中突变过的基因（突变频率 > 0）；
     ```python
     mutationFrequency = (mutMatrix.values != 0).sum(axis=1)
     mutMatrix = mutMatrix.loc[mutationFrequency > 0]
     ```

2. ✅ 是否对突变矩阵进行 gene 对齐不同：
   - main_adjusted.py 会将 mutMatrix 限制为与 expMatrix 的共同基因；
     ```python
     mutMatrix = mutMatrix.loc[commonGenes]
     ```
   - main_modified.py 注释掉此行，保留所有突变基因；
     ```python
     # mutMatrix = mutMatrix.loc[commonGenes]
     ```

3. ✅ 结果输出路径不同：
   - main_adjusted.py 输出到 results/PITCH_ad/<cancer>/genes_ranking.txt
   - main_modified.py 输出到 results/PITCH_md/<cancer>/genes_ranking.txt

4. ✅ （轻微差异）注释风格不同：
   - main_modified.py 注释中 `#` 后没有过多解释，代码注释略精简；
   - main_adjusted.py 注释更规范，有一定解释性描述。

⚠️ 其余部分（PPI 处理、z-score 归一化、构建超边、并行计算流程、调用 R 脚本等）在逻辑和实现上完全一致。
"""
import argparse
import os
import pandas as pd
import numpy as np
from scipy.stats import zscore
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

def preprocess_data(cancer):
    # Load data
    expMatrix = pd.read_csv(Path(__file__).resolve().parent.parent.parent/"data/processed_data"/cancer/"exp_tpm_data.tsv", sep='\t', index_col=0)
    mutMatrix = pd.read_csv(Path(__file__).resolve().parent.parent.parent/"data/processed_data"/cancer/"mut_data.tsv", sep='\t', index_col=0)
    PPI = pd.read_csv(Path(__file__).resolve().parent.parent.parent/'data/STRINGv12.txt', sep='\t')
    PPI.columns = ['Gene1', 'Gene2', 'Weight']
    Pathway = pd.read_csv('data/kegg_pathways_aggravated.txt', sep='\t', header=None)

    # Clean expression matrix gene names
    expMatrix.index = expMatrix.index.str.strip()
    expMatrix = expMatrix[~expMatrix.index.str.contains(r'[ .]')]
    expMatrix = expMatrix[expMatrix.index != '']

    # remain only columns with suffic "-01A"
    expMatrix = expMatrix.loc[:, expMatrix.columns.str.endswith('-01A')]
    mutMatrix = mutMatrix.loc[:, mutMatrix.columns.str.endswith('-01A')]
    # Filter mutation matrix
    numSamples = mutMatrix.shape[1]
    # mutationFrequency = (mutMatrix.values != 0).sum(axis=1) / numSamples
    # mutMatrix = mutMatrix.loc[mutationFrequency >= 0.02]
    mutationFrequency = (mutMatrix.values != 0).sum(axis=1)
    mutMatrix = mutMatrix.loc[mutationFrequency > 0]
    sampleMutationCounts = (mutMatrix.values != 0).sum(axis=0)
    mutMatrix = mutMatrix.loc[:, sampleMutationCounts >= 3]

    # Filter PPI by weight
    PPI = PPI[PPI['Weight'] >= 0.2]
    ppi_dict = build_ppi_dict(PPI)


    # Match samples and genes between expression and mutation matrices
    commonSamples = expMatrix.columns.intersection(mutMatrix.columns)
    mutMatrix = mutMatrix[commonSamples]
    expMatrix = expMatrix[commonSamples]
    mutMatrix.columns = expMatrix.columns = pd.Index(commonSamples).str.slice(0,12)
    commonGenes = expMatrix.index.intersection(mutMatrix.index)
    expMatrix = expMatrix.loc[commonGenes]
    # mutMatrix = mutMatrix.loc[commonGenes]

    # Z-score normalization by gene (row-wise)
    expMatrixZ = zscore_normalization(expMatrix, axis=1)

    # Identify differentially expressed genes
    diffGenes = {}
    allGenes = []
    threshold = 0.5
    for i in range(expMatrixZ.shape[1]):
        sample_i = expMatrixZ.columns[i]
        sampleZ = expMatrixZ.loc[:, sample_i]
        diffGenes[sample_i] = sampleZ.index[np.abs(sampleZ) > threshold].tolist() 
        allGenes.extend(diffGenes[sample_i])

    # allGenes = [gene for genes in diffGenes for gene in genes]
    uniqueGenes = list(set(allGenes))
    numUniqueGenes = len(uniqueGenes)
    print(f"Number of unique genes: {numUniqueGenes}")

    return expMatrix, mutMatrix, ppi_dict, Pathway, diffGenes

def build_ppi_dict(PPI):
    ppi_dict = {}
    for _, row in PPI.iterrows():
        g1, g2, w = row['Gene1'], row['Gene2'], row['Weight']
        key = tuple(sorted((g1, g2)))
        ppi_dict[key] = w

    return ppi_dict


def zscore_normalization(df, axis=0):
    """
    Perform z-score normalization on a pandas DataFrame.
    
    Parameters:
        df (pd.DataFrame): Input data.
        axis (int): 0 = z-score by column, 1 = z-score by row.
        
    Returns:
        pd.DataFrame: z-score normalized DataFrame.
    """
    if axis == 0:
        mean = df.mean(axis=0)
        std = df.std(axis=0)
        normalized = (df - mean) / std
    elif axis == 1:
        mean = df.mean(axis=1)
        std = df.std(axis=1)
        normalized = (df.sub(mean, axis=0)).div(std, axis=0)
    else:
        raise ValueError("Axis must be 0 (columns) or 1 (rows)")
    
    return normalized

def construct_hyperedges(allGenes, Pathway):
    hyperedges = []
    for _, row in Pathway.iterrows():
        if not isinstance(row[3], str):
            continue
        pathwayGenes = row[3].split(',')
        # mutGenes intersect sampleDEGs
        commonGenes = list(set(pathwayGenes).intersection(allGenes))
        if commonGenes:
            hyperedges.append(commonGenes)
    return hyperedges

def build_gene2hyperedges_map(hyperedges):
    gene2hyperedges = defaultdict(list)
    for idx, edge in enumerate(hyperedges):
        for gene in edge:
            gene2hyperedges[gene].append(idx)
    return gene2hyperedges

def compute_overlap_matrix(hyperedges):
    numHyperedges = len(hyperedges)
    overlapMatrix = np.zeros((numHyperedges, numHyperedges))

    for i in range(numHyperedges):
        for j in range(i + 1, numHyperedges):
            inter = max(0.1, len(set(hyperedges[i]).intersection(hyperedges[j])))
            min_len = min(len(hyperedges[i]), len(hyperedges[j]))
            overlapMatrix[i, j] = inter / min_len
            overlapMatrix[j, i] = overlapMatrix[i, j]

    return overlapMatrix

def compute_influence_for_sample(sampleIdx, mutMatrix, diffGenes, Pathway, ppi_dict):
    sampleName = mutMatrix.columns[sampleIdx]
    mutGenes = mutMatrix.index[mutMatrix.iloc[:, sampleIdx] == 1].tolist()
    sampleDEGs = diffGenes[sampleName]
    allGenes = list(set(mutGenes + sampleDEGs))
    hyperedges = construct_hyperedges(allGenes, Pathway)
    gene2hyperedges = build_gene2hyperedges_map(hyperedges)
    overlap_matrix = compute_overlap_matrix(hyperedges)

    numGenes = len(mutGenes)
    Infl = np.zeros(numGenes)

    for i, g_i in enumerate(mutGenes):
        for j, g_j in enumerate(allGenes):
            if i == j:
                continue
            # w_ij = PPI.loc[
            #     ((PPI['Gene1'] == g_i) & (PPI['Gene2'] == g_j)) |
            #     ((PPI['Gene1'] == g_j) & (PPI['Gene2'] == g_i)), 'Weight']
            # w_ij = w_ij.iloc[0] if not w_ij.empty else 0.01
            w_ij = ppi_dict.get(tuple(sorted((g_i, g_j))), 0.01)


            D_sum = 0
            for ep in gene2hyperedges.get(g_i, []):
                for eq in gene2hyperedges.get(g_j, []):
                    # inter = max(0.1,len(set(hyperedges[ep]).intersection(hyperedges[eq])))
                    # min_len = min(len(hyperedges[ep]), len(hyperedges[eq]))
                    D_sum += overlap_matrix[ep, eq]
            Infl[i] += D_sum * w_ij

    sorted_idx = np.argsort(-Infl)
    rankedGenes = [(mutGenes[idx], Infl[idx]) for idx in sorted_idx]
    return sampleName, rankedGenes

def run_driver_prioritization(cancer):
    expMatrix, mutMatrix, ppi_dict, Pathway, diffGenes = preprocess_data(cancer)

    args_list = [(i, mutMatrix, diffGenes, Pathway, ppi_dict) for i in range(mutMatrix.shape[1])]

    rankForSamples = {}
    with ProcessPoolExecutor(max_workers=61) as executor:
        futures = [executor.submit(compute_influence_for_sample, *args) for args in args_list]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Ranking samples"):
            sampleName, rankedGenes = future.result()
            rankForSamples[sampleName] = rankedGenes

    return rankForSamples

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Driver Gene Prioritization")
    parser.add_argument('-c',"--cancer", type=str, default='BRCA', help="Cancer type")
    args = parser.parse_args()
    cancer = args.cancer
    print(f"Running driver gene prioritization for {cancer}")
    
    all_ranks = run_driver_prioritization(cancer)

    results_dir = Path(__file__).resolve().parent.parent.parent/"results/PITCH_md"/cancer
    results_dir.mkdir(parents=True, exist_ok=True)
    result_file = results_dir / "genes_ranking.txt"
    with open(result_file, 'w') as f:
        for sample, ranked in all_ranks.items():
            gene_list = ','.join([gene for gene, score in ranked])
            f.write(f"{sample}\t{gene_list}\n")

    print(f"Finished processing {len(all_ranks)} samples.")
    
    pandas2ri.activate()
    
    r_source = robjects.r['source']
    r_source("../pernalized2cohort.R")
    personalized2cohort_score = robjects.r["personalized2cohort_score"]
    
    personalized2cohort_score(str(results_dir)+"/")