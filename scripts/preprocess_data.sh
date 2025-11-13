#!/bin/bash
# ==============================================================
# Script: preprocess_data.sh
# Purpose: Run full preprocessing pipeline for HyperNetWalk + WeSME
# Author: Xqxu
# Date: $(date +%Y-%m-%d)
# ==============================================================

set -e
set -o pipefail

# 初始化 conda 环境（关键行，允许脚本中使用 conda activate）
eval "$(conda shell.bash hook)"

# -----------------------------
# 1. 数据预处理（R 环境）
# -----------------------------
echo ">>> 激活 conda 环境: hypernetwalk"
conda activate hypernetwalk || { echo "无法激活 hypernetwalk 环境"; exit 1; }

echo ">>> 运行 R 数据预处理..."
Rscript src/data_preprocess.R
Rscript src/network_preprocess.R
echo ">>> 数据预处理完成。"

conda deactivate


# -----------------------------
# 2. WeSME 分析（Python 环境）
# -----------------------------
echo ">>> 进入 WeSME 工作目录"
cd src/wesme || { echo "目录 src/wesme 不存在"; exit 1; }

echo ">>> 激活 conda 环境: wesme"
conda activate wesme || { echo "无法激活 wesme 环境"; exit 1; }

# 定义癌症数据集数组
cancer_types=("PANCAN" "BRCA" "COAD" "HNSC" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "PRAD" "STAD" "THCA" "UCEC")
# cancer_types=("BRCA" "COAD" "HNSC" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC")
# -----------------------------
# 2.1 预处理突变矩阵并计算 ME 值
# -----------------------------
echo ">>> 开始并行运行 WeSME 主分析..."
printf "%s\n" "${cancer_types[@]}" | xargs -P 8 -I {} bash -c '
    echo "--- Processing cancer type: {} ---"
    python mut_mat_to_list.py {} smut &&
    python comp_sample_weights.py {} smut &&
    python run_weighted_sampling.py {} smut 1000 &&
    python comp_me_for_all_pairs.py {} smut 0.1
'

# -----------------------------
# 2.2 运行置换检验
# -----------------------------
echo ">>> 开始并行运行置换检验..."
for cancer in "${cancer_types[@]}"; do
    echo "--- Running permutations for $cancer ---"
    seq 0 99 | xargs -P 100 -I {} bash -c "
        cd . && \
        echo 当前目录: \$(pwd)
        ls preproc/mut/$cancer/
        python run_permute_data.py $cancer smut {} 100 &&
        python comp_me_for_all_pairs.py $cancer smut 0.1 \
            -m preproc/permuted_mut/$cancer/${cancer}_smut_permuted_cover_{}.txt \
            -o preproc/permuted_pv/$cancer/${cancer}_smut_perm_pv_{}.txt
    "
done

# -----------------------------
# 2.3 计算 FDR
# -----------------------------
echo ">>> 计算 FDR..."
printf "%s\n" "${cancer_types[@]}" | xargs -P 8 -I {} bash -c '
    python comp_fdr.py {} smut me 100 2
'

# -----------------------------
# 3. 收尾
# -----------------------------
conda deactivate
echo ">>> 已返回 hypernetwalk 环境"
conda activate hypernetwalk || echo "警告：hypernetwalk 环境未重新激活"

echo "✅ 所有步骤执行完毕。"
