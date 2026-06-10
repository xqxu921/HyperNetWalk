#!/bin/bash

# 16 cancer types
# cancers=("BLCA" "COAD" "ESCA" "GBM" "HNSC" "KIRC" "LGG" 
#          "LIHC" "LUAD" "LUSC" "PAAD" "SKCM" "STAD" "THCA" "UCEC")
# cancers=("KIRP" "PRAD" "READ")
# cancers=("PANCAN")
cancers=("BRCA" "COAD" "HNSC" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "PRAD" "STAD" "THCA" "UCEC")
for cancer in "${cancers[@]}"
do
    echo "Processing $cancer..."

    # python main.py -c "$cancer"
    python main_adjusted.py -c "$cancer"

    echo "$cancer done."
done

echo "All cancers processed!"