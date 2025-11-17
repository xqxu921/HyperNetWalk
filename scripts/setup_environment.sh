#!/bin/bash
set -e

echo "=== 初始化 conda 环境 ==="

# 1. 尝试自动检测 conda（适用于所有用户）
if command -v conda >/dev/null 2>&1; then
    echo "检测到 conda：$(command -v conda)"
else
    # 2. 尝试从常见安装路径加载 conda
    for p in \
        "$HOME/miniconda3/etc/profile.d/conda.sh" \
        "$HOME/anaconda3/etc/profile.d/conda.sh" \
        "/opt/miniconda3/etc/profile.d/conda.sh" \
        "/opt/anaconda3/etc/profile.d/conda.sh"
    do
        if [ -f "$p" ]; then
            echo "从 $p 加载 conda"
            source "$p"
            break
        fi
    done
fi

# 3. 再次检查 conda 是否可用
if ! command -v conda >/dev/null 2>&1; then
    echo "错误：无法找到 conda，请安装 Miniconda 或 Anaconda 后再运行此脚本。"
    exit 1
fi

echo "=== conda 初始化完成 ==="

echo "=== [1/3] 创建 Python 环境 (WeSME) ==="
conda env create -f src/wesme/environment.yml || conda env update -f src/wesme/environment.yml
echo "Python 环境 wesme 创建完成。"

echo "=== [2/3] 创建 R 环境 (HyperNetWalk) ==="
conda env create -f environment.yml || (conda env remove -n hypernetwalk && conda env create -f environment.yml)
echo "R 环境 hypernetwalk 创建完成。"

echo "=== [3/3] 初始化 R 包依赖 (renv) ==="
conda run -n hypernetwalk Rscript -e 'renv::restore()'
echo "R 包依赖恢复完成。"

echo "=== 环境已准备就绪！ ==="