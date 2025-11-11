#!/bin/bash
set -e

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