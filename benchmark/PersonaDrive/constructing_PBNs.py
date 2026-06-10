#***********************************
# Run the command 'python constructing_personalized_networks.py' in the Terminal to start the
# the construction of the personalized networkxs
# Or directly run this file in your chosen IDE.
# This may take several hours.
# When is done, the files of prepared data are saved in the subfolders of graphs/[dataset]/[cancer].
#***********************************
import argparse
from tqdm import tqdm as tqdm
import os
import sys
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
from pathlib import Path

## Load DEGs matrix
def load_DEGs(file):
    DEGs_matrix = pd.read_csv(file, index_col=0,sep=",")
    DEGs_matrix = DEGs_matrix.transpose()
    # DEGs_matrix=DEGs_matrix.rename(columns = {i:i.replace('.','-') for i in DEGs_matrix},index = {i:i.replace('.','-') for i in DEGs_matrix.index})
    return DEGs_matrix

## Load MUTATION matrix
def load_mutations(file):
    M_matrix = pd.read_csv(file, index_col=0, sep="\t")
    M_matrix = M_matrix.rename(columns = {i:i.replace('.','-') for i in M_matrix})
    return M_matrix

## Load PPI
def load_ppi(file):
    ppi_matrix = pd.read_csv(file,sep="\t")
    ppi = {}
    columns_name = ppi_matrix.columns
    ppi_matrix = ppi_matrix[ppi_matrix[columns_name[2]]>=0.7].reset_index(drop=True)
    for i in range(len(ppi_matrix.index)):
        ge1 = ppi_matrix.loc[i,columns_name[0]]
        ge2 = ppi_matrix.loc[i,columns_name[1]]
        if ge1 not in ppi:
            ppi[ge1]=set()
            ppi[ge1].add(ge2)
        else:
            ppi[ge1].add(ge2)
        if ge2 not in ppi:
            ppi[ge2]=set()
            ppi[ge2].add(ge1)
        else:
            ppi[ge2].add(ge1)
    return ppi

## construct personalized networks
def construct_PBNs(M,O,PPI,cancer,dataset,network):

    #create the output folder
    if not os.path.exists("graphs/"+dataset+"/"):
        os.mkdir("graphs/"+dataset+"/")
    if not os.path.exists("graphs/"+dataset+"/"+cancer+"_"+network):
        os.mkdir("graphs/"+dataset+"/"+cancer+"_"+network)

    #####################################
    #   Starting with mutation data     #
    #####################################
    #list samples in M matrix
    samples = M.columns.tolist()
    M = M.isin([1])

    #dictionary that contains samples ids and their mutations
    M_dictionary ={}
    for i in samples:
        m_i = M.index[M[i]].tolist()
        M_dictionary[i]= list(set(m_i))

    #####################################
    #   Pre-processing outliers data    #
    #####################################

    #genes is a set of genes in outliers dataset
    genes = O.index.tolist()
    #list of outliers
    D_set= []
    for i in tqdm(samples):
        for gene in genes:
            if i in O.columns:
                if O[i][gene] == True:
                    # outliers are in form of: CCLEXXXX_TP53
                    D_set.append(str(i+'_'+gene))
    print("There are: ",len(D_set)," outliers.")



    #####################################
    #  Constructing Bipartite networks  #
    #####################################

    for i in tqdm(samples):

        if 'BP_'+i+'.gml' not in os.listdir("graphs/"+dataset+"/"+cancer+"_"+network):

            # initialize the bipartite networks B_i
            Bi = nx.Graph()
            Bi.add_nodes_from(M_dictionary[i], bipartite=0)
            Bi.add_nodes_from(D_set, bipartite=1)


            Mi=M_dictionary[i]
            for u in tqdm(Mi):
                if u in PPI:
                    for v in D_set:
                        v_d = v.split('_')[1]
                        # check if the selected mutated gene is mutated in the sample of current outlier
                        if u in M_dictionary[v.split('_')[0]]:
                            if v_d in PPI:
                                #check if outlier is a neighbor of the mutation
                                if v_d in PPI[u]:
                                    if v_d !=u:
                                        Bi.add_edge(u,v)


            bipartite0_nodes = {n for n, d in Bi.nodes(data=True) if d['bipartite']==0}
            bipartite1_nodes = {n for n, d in Bi.nodes(data=True) if d['bipartite']==1}

            #deleting nodes with 0 degree
            for v in bipartite1_nodes:
                if int(Bi.degree(v))==0:
                    Bi.remove_node(v)

            for u in bipartite0_nodes:
                if int(Bi.degree(u))==0:
                    Bi.remove_node(u)


            bipartite0_nodes = {n for n, d in Bi.nodes(data=True) if d['bipartite']==0}
            bipartite1_nodes = {n for n, d in Bi.nodes(data=True) if d['bipartite']==1}

            print("Mutations:  ",len(bipartite0_nodes),"   Outliers:",len(bipartite1_nodes))
            nx.write_gml(Bi, "graphs/"+dataset+"/"+cancer +"_"+network + '/BP_'+i+'.gml')


if __name__ == "__main__":
    data_folder="data/"

    description="Personalized Bipartite Networks (PBNs)"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-d', "--dataset", type=str, required=False, default='TCGA', help="Dataset: TCGA or CCLE")
    parser.add_argument('-c', '--cancer', type=str, required=False, default="COAD", help="Cancer Type")
    parser.add_argument('-n', '--network', type=str, required=False, default="ST11", help="PPI Network")

    args = parser.parse_args()
    dataset = args.dataset
    cancer = args.cancer
    network = args.network
    #~~~~~~~~~~~~~Step 1：~~~~~~~~~~~~~~~~~~
    # load PPI data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    networks={
        "ST11":"data/string_v11.5_network.csv",
        "ST10":"data/string_v10.5_network.csv",
        "DW":"data/dawnrank_network.csv",
        "ST12":Path(__file__).resolve().parent.parent.parent / "data/STRINGv12.txt"
    }
    ppi_file = networks[network]
    ppi = load_ppi(ppi_file)
    print('PPI network loaded...')

    #~~~~~~~~~~~~~Step 2：~~~~~~~~~~~~~~~~~~
    # load outliers data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # DEGs_matrix = load_DEGs(data_folder+cancer+"_"+dataset+"/DEGs.csv")
    DEGs_matrix = load_DEGs(Path(__file__).resolve().parent.parent.parent / "data/processed_data" / cancer / "DEGs.csv")
    print('DEGs loaded...')

    #~~~~~~~~~~~~~Step 3：~~~~~~~~~~~~~~~~~~
    # load mutation data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # M_matrix = load_mutations(data_folder+cancer+"_"+dataset+"/MUT.csv")
    M_matrix = load_mutations(Path(__file__).resolve().parent.parent.parent / "data/processed_data" / cancer / "mut_data.tsv")
    print('Mutations loaded...')

    #~~~~~~~~~~~~~Step 4：~~~~~~~~~~~~~~~~~~
    # pre-process the matrices to contain only genes in the PPI
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEGs_matrix = DEGs_matrix.drop(index=[n for n in DEGs_matrix.index if n not in ppi])
    M_matrix = M_matrix.drop(index=[n for n in M_matrix.index if n not in ppi])
    
    M_matrix = M_matrix.loc[:, M_matrix.columns.str.endswith(("-01A", "-01"))]
    common_samples = list(set(DEGs_matrix.columns).intersection(set(M_matrix.columns)))
    DEGs_matrix = DEGs_matrix[common_samples]
    M_matrix = M_matrix[common_samples]
    M_matrix.columns = DEGs_matrix.columns = pd.Index(common_samples).str.slice(0, 12)

    #~~~~~~~~~~~~~Step 5：~~~~~~~~~~~~~~~~~~
    # construct the personalized networks
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    construct_PBNs(M_matrix,DEGs_matrix,ppi,cancer,dataset,network)