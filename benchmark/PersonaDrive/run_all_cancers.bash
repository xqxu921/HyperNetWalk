#!/bin/bash

# 16 cancer types
cancers=("BRCA" "COAD" "HNSC" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "PRAD" "STAD" "THCA" "UCEC")
# cancers=("KIRP" "PRAD" "READ")
# cancers=("PANCAN")
for cancer in "${cancers[@]}"
do
    echo "Processing $cancer..."

    # python constructing_PBNs.py -d TCGA -c "$cancer" -n ST12
    python PersonaDrive.py -d TCGA -c "$cancer" -n ST12

    echo "$cancer done."
done

echo "All cancers processed!"