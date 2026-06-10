#***********************************
# Run the command 'python PersonaDrive.py' in the Terminal.
# Or directly run this file in the IDE.

# Note that, before running this script, you should run the script of constructing_personalized_networks.py to prepare the data for PersonaDrive.
#
#

# This may take several minutes to finish all the samples.
# Exact time mainly depends on your computer and the number of samples in your dataset.
#***********************************
import argparse
from tqdm import tqdm as tqdm
import os
import sys
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
import math
from pathlib import Path
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter

#Load pathways
def load_pathways():
    with open("data/kegg_pathways_v1.txt", 'r')as ifile:
        pathways={}
        for line in ifile:
            line = line.rstrip('\n')
            p=line.split('\t\t\t')[0] #pathway ID
            genes=set(line.split('\t\t\t')[1].split(",") [:-1])
            for g in genes:
                if g not in pathways:
                    pathways[g]=set()
                    pathways[g].add(p)
                else:
                    pathways[g].add(p)
    return pathways

def calculate_pps(graphs,dataset,cancer,network):
    samples=[sample.split('_')[-1].replace('.gml','') for sample in os.listdir(graphs+dataset+"/"+cancer+"_"+network) if '.gml' in sample]

    #load graphs
    D_data={}

    for i in tqdm(os.listdir(graphs+dataset+"/"+cancer+"_"+network)): # for each sample i in S
        if '.gml' in i:
            sample = i.split('_')[-1].replace('.gml','') # retreive samle ID from the graph name

            Bi = nx.read_gml(graphs+dataset+"/"+cancer+"_"+network+"/"+i) #load Bi graph
            # set of all outliers O_final in the second partition connected to M_i
            D_Bi = [n for n, d in Bi.nodes(data=True) if d['bipartite']==1]
            Di=[]

            # select onmy outliers that corresponds to samlpe_i
            for v_d in D_Bi:
                if v_d.split("_")[0]==sample:
                    Di.append(v_d.split("_")[1])

            D_data[sample]=set(Di)
    pps=pd.DataFrame(0.0,columns=samples,index=samples)
    for i in pps: #for each sample i in S
        D_i=D_data[i]
        for j in pps: #for each sample i in S
            D_j=D_data[j]
            #calculate the similarity score
            if len(D_i)*len(D_j)!=0:
                pps.loc[j,i]=math.pow(len(D_i.intersection(D_j)),2)/(len(D_i)*len(D_j))
    pps.to_csv(Path(__file__).resolve().parent.parent.parent / "results/PersonaDrive" / cancer / "pps_matrix.csv")

# calculate influence scores
def calculate_infl_scores(Bi,Mi,sample_i,similarity,pathways):
    influence_scores={}
    for u in Mi:
        influence_scores[u]=0.0
        for v in Bi.neighbors(u):
            sample_j=v.split("_")[0]
            # check that the mutaated gene u and the outlier gene v exist in the pathways data
            if v.split("_")[1] in pathways and u in pathways:
                # calculate cp score
                ppc_u_v=len(set(pathways[v.split("_")[1]]).intersection(set(pathways[u])))
                # retreive the patient similarity score
                pps_i_j=similarity.loc[sample_j,sample_i]

                influence_scores[u]+= ppc_u_v * pps_i_j

    return influence_scores

#prioritize mutated genes
def PersonaDrive(graphs,pps_matrix,cancer,dataset,pathways,network):
    # with open("results/"+dataset+"/"+cancer+"_"+network+"/PersonaDrive.txt","w") as ofile:
    with open(Path(__file__).resolve().parent.parent.parent / "results/PersonaDrive_hrw" / cancer /"genes_ranking.txt","w") as ofile:    
        personalized_drivers={}
        scores = {}
        for i in tqdm(os.listdir(graphs+dataset+"/"+cancer+"_"+network)):
            if '.gml' in i:
                #get sample id
                sample = i.split('_')[-1].replace('.gml','')

                #load the graph file
                Bi = nx.read_gml(graphs+dataset+"/"+cancer+"_"+network+"/"+i)
                #list all mutations
                Mi = [n for n, d in Bi.nodes(data=True) if d['bipartite']==0]
                if len(Mi) == 0:
                    continue
                else:
                    #assign edge weights
                    infl_scores = calculate_infl_scores(Bi, Mi,sample,pps_matrix,pathways)
                    drivers=sorted(infl_scores.items(), key=lambda kv: (-kv[1], kv[0]))
                    ofile.write(sample+"\t")
                    # for u in [u_i[0] for u_i in drivers]:
                    #     ofile.write(u+"\t")
                    # ofile.write("\n")
                    ofile.write(",".join([u_i[0] for u_i in drivers]) + "\n")
                    scores[sample]=drivers
    return scores


if __name__ == "__main__":
    graphs="graphs/"

    description="Ranking Mutated Genes in PBNs"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-d', "--dataset", type=str, required=False, default='TCGA', help="Dataset: TCGA or CCLE")
    parser.add_argument('-c', '--cancer', type=str, required=False, default="COAD", help="Cancer Type")
    parser.add_argument('-n', '--network', type=str, required=False, default="ST11", help="PPI Network")

    args = parser.parse_args()
    dataset = args.dataset
    cancer = args.cancer
    network = args.network


    #~~~~~~~~~~~~~Step 1：~~~~~~~~~~~~~~~~~~
    # load pathways data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pathways = load_pathways()
    print('pathways loaded...')


    #~~~~~~~~~~~~~Step 1：~~~~~~~~~~~~~~~~~~
    # construct pps matrix
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    try:
        # os.mkdir("results/"+dataset+"/"+cancer+"_"+network)
        if not os.path.exists(Path(__file__).resolve().parent.parent.parent / "results/PersonaDrive_hrw"):
            os.mkdir(Path(__file__).resolve().parent.parent.parent / "results/PersonaDrive_hrw")
        os.mkdir(Path(__file__).resolve().parent.parent.parent / "results/PersonaDrive_hrw" / cancer)
    except:
        print("Folder Exist!")
    if "pps_matrix.csv" not in os.listdir(Path(__file__).resolve().parent.parent.parent / "results/PersonaDrive" / cancer):
        calculate_pps(graphs,dataset,cancer,network)
    pps_file=Path(__file__).resolve().parent.parent.parent / "results/PersonaDrive" / cancer / "pps_matrix.csv"
    pps_matrix=pd.read_csv(pps_file,index_col=0)
    print('PPS matrix loaded...')

    #~~~~~~~~~~~~~Step 1：~~~~~~~~~~~~~~~~~~
    # run PersonaDrive
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scores = PersonaDrive(graphs,pps_matrix,cancer,dataset,pathways,network)
    #把scores这个dict转为矩阵，注意每个元素代表一个样本，每个样本对应不同长度的基因打分tuples，对应矩阵的行表示所有基因并集，列表示样本，元素值为对应基因在对应样本中的打分，如果该基因不在该样本中，则为0
    all_genes = set()
    for sample in scores:
        for gene, score in scores[sample]:
            all_genes.add(gene)
    all_genes = sorted(list(all_genes))
    gene_score_matrix = pd.DataFrame(0.0, index=all_genes, columns=scores.keys())
    for sample in scores:
        for gene, score in scores[sample]:
            gene_score_matrix.loc[gene, sample] = score

    # save the gene score matrix
    gene_score_matrix.to_csv(Path(__file__).resolve().parent.parent.parent / "results/PersonaDrive" / cancer / "gene_score_matrix.csv")
            
    # Use a local conversion context instead of the deprecated pandas2ri.activate()
    r_source = ro.r['source']
    r_source("../pernalized2cohort.R")
    personalized2cohort_score_new = ro.r["personalized2cohort_score_new"]
    #把gene_score_matrix转为R的Matrix格式，行名为基因，列名为样本，元素值为对应的打分
    with localconverter(default_converter + pandas2ri.converter):
        r_gene_score_matrix = pandas2ri.py2rpy(gene_score_matrix)
    # Call the R function inside a localconverter so pandas<->R data.frame
    # conversions are applied only within this block (no global activation).
    with localconverter(default_converter + pandas2ri.converter):
        _ = personalized2cohort_score_new(r_gene_score_matrix, str(Path(__file__).resolve().parent.parent.parent / "results/PersonaDrive_hrw" / cancer / ""))


    # personalized2cohort_score = ro.r["personalized2cohort_score"]

    # input_file = Path(__file__).resolve().parent.parent.parent / "results/PersonaDrive" / cancer
    # # Call the R function inside a localconverter so pandas<->R data.frame
    # # conversions are applied only within this block (no global activation).
    # with localconverter(default_converter + pandas2ri.converter):
    #     _ = personalized2cohort_score(str(input_file) + "/")