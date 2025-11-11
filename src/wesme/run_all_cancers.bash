!/bin/bash
定义癌症数据集数组
cancer_types=("BRCA" "BLCA" "COAD" "ESCA" "GBM" "HNSC" "KIRC" "KIRP" "LGG"
              "LIHC" "LUAD" "LUSC" "PAAD" "PRAD" "READ" "SKCM" "STAD" "THCA" "UCEC")
cancer_types=("BRCA" "COAD" "HNSC" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "PRAD" "STAD" "THCA" "UCEC")
# 使用xargs并行运行，-P指定并行进程数
# printf "%s\n" "${cancer_types[@]}" | xargs -P 8 -I {} bash -c "
#     python mut_mat_to_list.py {} smut &&
#     python comp_sample_weights.py {} smut &&
#     python run_weighted_sampling.py {} smut 1000 &&
#     python comp_me_for_all_pairs.py {} smut 0.1
# "

for cancer in BRCA COAD HNSC KIRC KIRP LIHC LUAD LUSC PRAD STAD THCA UCEC; do
    seq 0 99 | xargs -P 100 -I {} bash -c "
        python run_permute_data.py $cancer smut {} 100 &&
        python comp_me_for_all_pairs.py $cancer smut 0.1 \
        -m preproc/permuted_mut/$cancer/${cancer}_smut_permuted_cover_{}.txt \
        -o preproc/permuted_pv/$cancer/${cancer}_smut_perm_pv_{}.txt
    "
done

printf "%s\n" "${cancer_types[@]}" | xargs -P 8 -I {} bash -c "
    python comp_fdr.py {} smut me 100 2
"