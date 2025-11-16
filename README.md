# 🔬 HyperNetWalk

**HyperNetWalk** is a hypergraph-based framework for **pan-cancer driver gene identification** across multi-omics layers.  
It integrates **protein–protein interaction (PPI)**, **gene regulatory (GRN)**, and **mutual exclusivity (ME)** networks to identify driver genes at both **cohort** and **individual** levels.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R >= 4.5.2](https://img.shields.io/badge/R-%3E%3D4.5.2-green.svg)](https://cran.r-project.org/)
[![Conda](https://img.shields.io/badge/environment-conda-orange.svg)](https://docs.conda.io/)
[![WeSME](https://img.shields.io/badge/Python-WeSME-yellow.svg)](https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/index.cgi#wesme)

---

## 📋 Table of Contents

- [1. Environment Setup](#1-environment-setup)
- [2. Data Preparation](#2-data-preparation)
- [3. Running HyperNetWalk](#3-running-hypernetwalk)
- [4. Results Evaluation](#4-results-evaluation)
- [5. Quick Validation (Recommended)](#5-quick-validation-recommended)
- [6. Repository Structure](#6-repository-structure)
- [7. Citation](#7-citation)
- [8. 中文说明(简要)](#8-中文说明简要)

---

## 🧩 1. Environment Setup

Clone this repository and automatically configure all environments:

```bash
git clone https://github.com/xqxu921/HyperNetWalk.git
cd HyperNetWalk
bash scripts/setup_environment.sh
```

**This script will:**
- ✅ Create the **Python environment** (for WeSME preprocessing)
- ✅ Create the **R environment** (`hypernetwalk`, R 4.5.2 with proper `crossprod()` behavior)
- ✅ Restore all R dependencies using `renv::restore()`

After installation, activate the conda environment:

```bash
conda activate hypernetwalk
```

---

## 📦 2. Data Preparation

You can choose **one of the following methods**:

### Option 1: Download Preprocessed Data from SourceForge (Recommended)

#### Method 1: Download ZIP Files (Recommended) ⭐

Download the following ZIP files from SourceForge:
- `DRIVER.zip` - Driver gene annotations
- `NETWORK.zip` - PPI, GRN networks and mutual exclusive networks
- `metadata.zip` - Gene metadata
- `processed.zip` - Preprocessed omics data

**Using wget:**

```bash
# Download all ZIP files
wget https://sourceforge.net/projects/hypernetwalk/files/data/DRIVER.zip/download -O DRIVER.zip
wget https://sourceforge.net/projects/hypernetwalk/files/data/NETWORK.zip/download -O NETWORK.zip
wget https://sourceforge.net/projects/hypernetwalk/files/data/metadata.zip/download -O metadata.zip
wget https://sourceforge.net/projects/hypernetwalk/files/data/processed.zip/download -O processed.zip

# Extract all files to data directory
unzip -q DRIVER.zip -d data/
unzip -q NETWORK.zip -d data/
unzip -q metadata.zip -d data/
unzip -q processed.zip -d data/

# Clean up ZIP files (optional)
rm DRIVER.zip NETWORK.zip metadata.zip processed.zip
```

**Using curl:**

```bash
# Download all ZIP files
curl -L https://sourceforge.net/projects/hypernetwalk/files/data/DRIVER.zip/download -o DRIVER.zip
curl -L https://sourceforge.net/projects/hypernetwalk/files/data/NETWORK.zip/download -o NETWORK.zip
curl -L https://sourceforge.net/projects/hypernetwalk/files/data/metadata.zip/download -o metadata.zip
curl -L https://sourceforge.net/projects/hypernetwalk/files/data/processed.zip/download -o processed.zip

# Extract all files to data directory
unzip -q DRIVER.zip -d data/
unzip -q NETWORK.zip -d data/
unzip -q metadata.zip -d data/
unzip -q processed.zip -d data/

# Clean up ZIP files (optional)
rm DRIVER.zip NETWORK.zip metadata.zip processed.zip
```

#### Method 2: Manual Download

1. Visit: https://sourceforge.net/projects/hypernetwalk/files/data/
2. Click to download each ZIP file:
   - `DRIVER.zip`
   - `NETWORK.zip`
   - `metadata.zip`
   - `processed.zip`
3. Extract all ZIP files to the `data/` directory in your project

**Expected directory structure after extraction:**
```
data/
├── DRIVER/
│   ├── CGC_Tier1.tsv
│   └── Compendium_Cancer_Genes.tsv
├── NETWORK/
│   ├── STRINGv12.txt
│   ├── RegNet_human_V2.txt
│   ├── BRCA_me_net.txt
│   └── [ME net of other cancer types]
├── metadata/
│   ├── gencode.v36.annotation.gtf.gene.probemap
│   └── 9606.protein.info.v12.0.txt
└── processed/
    ├── PANCAN/
    ├── BRCA/
    └── [other cancer types]
```
**⚠️ Important Notes:**
- Ensure you have `unzip` installed on your system
- Total download size: Check SourceForge for file sizes
- Ensure you have sufficient disk space for both ZIP files and extracted data

### Option 2: Download Raw Data and Preprocess Locally

```bash
bash scripts/download_raw_data.sh
bash scripts/preprocess_data.sh
```

**This will obtain:**
- 🧬 TCGA Pan-cancer expression and mutation raw files
- 🎗️ BRCA and 11 other cancer types multi-omics raw data
- 🔗 STRING v12 protein-protein interaction network raw files
- 📊 RegNetwork transcriptional regulatory network raw files
- 📚 COSMIC CGC Tier1 and IntOGen driver gene annotations raw files

---

## 🚀 3. Running HyperNetWalk

### Create Necessary Directories

```bash
mkdir -p results logs
```

### Manual Step-by-Step Execution

#### (1) Pan-cancer Cohort-level Prediction

```bash
/usr/bin/time -v -o logs/pancan_resource_usage.txt \
  Rscript src/run_hypernetwalk.R \
    --mode pancancer \
    --level cohort \
    --input data/processed \
    --ppi data/NETWORK/STRINGv12.txt \
    --grn data/NETWORK/RegNet_human_V2.txt \
    --output results/ \
    --cores 64
```

#### (2) Single Cancer Cohort-level Prediction (e.g., BRCA)

```bash
/usr/bin/time -v -o logs/brca_cohort_resource.txt \
  Rscript src/run_hypernetwalk.R \
    --mode single_cancer \
    --level cohort \
    --cancer_type BRCA \
    --input data/processed/ \
    --ppi data/NETWORK/STRINGv12.txt \
    --grn data/NETWORK/RegNet_human_V2.txt \
    --output results/ \
    --cores 64
```

#### (3) Single Cancer Individual-level Prediction (e.g., BRCA)

```bash
/usr/bin/time -v -o logs/brca_individual_resource.txt \
  Rscript src/run_hypernetwalk.R \
    --mode single_cancer \
    --level individual \
    --cancer_type BRCA \
    --input data/processed/ \
    --ppi data/NETWORK/STRINGv12.txt \
    --grn data/NETWORK/RegNet_human_V2.txt \
    --output results/ \
    --cores 64
```

---

## 📊 4. Results Evaluation

### Evaluate Pan-cancer Cohort Results

```bash
Rscript src/evaluation.R \
  --mode pancancer \
  --level cohort \
  --predicted results/PANCAN \
  --benchmark CGC \
  --output results/PANCAN/evaluation_results.txt
```

### Evaluate BRCA Cohort Results

```bash
Rscript src/evaluation.R \
  --mode single_cancer \
  --level cohort \
  --cancer_type BRCA \
  --predicted results/BRCA/ \
  --benchmark CGC \
  --output results/BRCA/evaluation_results.txt
```

### Evaluate BRCA Individual Results

```bash
Rscript src/evaluation.R \
  --mode single_cancer \
  --level individual \
  --cancer_type BRCA \
  --predicted results/BRCA/ \
  --benchmark CGC \
  --output results/BRCA/evaluation_results.txt
```

---

## ⚡ 5. Quick Validation (Recommended)

### Step 1: Ensure Proper Installation

Make sure you are in the HyperNetWalk directory:

```bash
cd /path/to/HyperNetWalk
```

Activate the conda environment:

```bash
conda activate hypernetwalk
```

Ensure scripts are executable:

```bash
chmod +x scripts/run_all_tests.sh scripts/evaluate_all_results.sh
```

### Step 2: Download Preprocessed Data from SourceForge

Download the required data and scripts from:  
👉 [https://sourceforge.net/projects/hypernetwalk/files/data/](https://sourceforge.net/projects/hypernetwalk/files/data/)

### Step 3: Run Complete Testing Workflow

```bash
bash scripts/run_all_tests.sh
```

**This script will sequentially run:**
- Pan-cancer cohort prediction
- 12 cancer types predictions (both cohort and individual levels for each cancer type)

### Step 4: View Summary Results

After completion, evaluate all results:

```bash
bash scripts/evaluate_all_results.sh
cat results/summary_report.txt
less results/detailed_report.txt
```

### Step 5: Check Resource Usage

View resource usage logs:

```bash
ls -lh logs/
cat logs/*_resource_usage.txt | grep "Maximum resident set size"
```

---

## 📁 6. Repository Structure

```
HyperNetWalk/
├── data/                          # Omics and network data
│   ├── DRIVER/                    # Driver gene annotations
│   ├── NETWORK/                   # PPI and GRN networks
│   │   ├── STRINGv12.txt         # STRING v12 PPI network
│   │   └── RegNet_human_V2.txt   # RegNetwork GRN
│   ├── metadata/                  # Sample metadata
│   └── processed/                 # Preprocessed omics data
├── results/                       # Model outputs
│   ├── PANCAN/                    # Pan-cancer results
│   ├── BRCA/                      # BRCA results
│   ├── summary_report.txt         # Summary report
│   └── detailed_report.txt        # Detailed report
├── logs/                          # Resource usage logs
├── scripts/                       # Setup & automation scripts
│   ├── setup_environment.sh       # Environment configuration
│   ├── download_raw_data.sh       # Data download script
│   ├── preprocess_data.sh         # Data preprocessing script
│   ├── run_all_tests.sh           # Automated testing workflow
│   └── evaluate_all_results.sh    # Automated evaluation workflow
├── src/                           # Source code
│   ├── wesme/                     # Python-based preprocessing module
│   ├── run_hypernetwalk.R         # Main HyperNetWalk model
│   └── evaluation.R               # Evaluation and benchmarking
├── environment.yml                # Conda environment for R
├── renv.lock                      # R package snapshot
└── README.md                      # This file
```

---

## ✨ 7. Citation

If you use **HyperNetWalk** in your research, please cite:

```bibtex
@article{xu2025hypernetwalk,
  title={HyperNetWalk: Integrative Hypergraph-based Framework for Pan-cancer Driver Gene Identification},
  author={Xu, XQ and others},
  journal={Journal Name},
  year={2025},
  publisher={Publisher}
}
```

---

## 🇨🇳 8. 中文说明(简要)

**HyperNetWalk** 是一个基于超图的泛癌驱动基因识别框架，整合了多组学数据层。

### 快速开始

#### 一、配置运行环境

```bash
# 克隆仓库并配置环境
git clone https://github.com/xqxu921/HyperNetWalk.git
cd HyperNetWalk
bash scripts/setup_environment.sh

# 激活conda环境
conda activate hypernetwalk
```

#### 二、数据准备

**方式一：从SourceForge下载预处理数据（推荐）**

访问：https://sourceforge.net/projects/hypernetwalk/files/data/

**使用 wget 下载（推荐）：**
```bash
# 下载所有ZIP文件
wget https://sourceforge.net/projects/hypernetwalk/files/data/DRIVER.zip/download -O DRIVER.zip
wget https://sourceforge.net/projects/hypernetwalk/files/data/NETWORK.zip/download -O NETWORK.zip
wget https://sourceforge.net/projects/hypernetwalk/files/data/metadata.zip/download -O metadata.zip
wget https://sourceforge.net/projects/hypernetwalk/files/data/processed.zip/download -O processed.zip

# 解压ZIP文件到data/目录下
unzip -q DRIVER.zip -d data/
unzip -q NETWORK.zip -d data/
unzip -q metadata.zip -d data/
unzip -q processed.zip -d data/

# 清除ZIP文件(可选)
rm DRIVER.zip NETWORK.zip metadata.zip processed.zip
```

**使用 curl：**
```bash
# 下载所有ZIP文件
curl -L https://sourceforge.net/projects/hypernetwalk/files/data/DRIVER.zip/download -o DRIVER.zip
curl -L https://sourceforge.net/projects/hypernetwalk/files/data/NETWORK.zip/download -o NETWORK.zip
curl -L https://sourceforge.net/projects/hypernetwalk/files/data/metadata.zip/download -o metadata.zip
curl -L https://sourceforge.net/projects/hypernetwalk/files/data/processed.zip/download -o processed.zip

# 解压ZIP文件到data/目录下
unzip -q DRIVER.zip -d data/
unzip -q NETWORK.zip -d data/
unzip -q metadata.zip -d data/
unzip -q processed.zip -d data/

# 清除ZIP文件(可选)
rm DRIVER.zip NETWORK.zip metadata.zip processed.zip
```

**手动下载：**
下载以下目录：`DRIVER.zip`, `NETWORK.zip`, `metadata.zip`, `processed.zip`，并解压放置到项目的 `data/` 目录下

**方式二：使用脚本下载原始数据并预处理**

```bash
bash scripts/download_raw_data.sh
bash scripts/preprocess_data.sh
```

#### 三、运行HyperNetWalk

```bash
# 创建必要的目录
mkdir -p results logs

# Pan-cancer群体水平预测
/usr/bin/time -v -o logs/pancan_resource_usage.txt \
  Rscript src/run_hypernetwalk.R \
    --mode pancancer \
    --level cohort \
    --input data/processed \
    --ppi data/NETWORK/STRINGv12.txt \
    --grn data/NETWORK/RegNet_human_V2.txt \
    --output results/ \
    --cores 64
```

#### 四、快速验证（推荐）

```bash
# 1. 确保脚本可执行
chmod +x scripts/run_all_tests.sh scripts/evaluate_all_results.sh

# 2. 运行完整测试流程
bash scripts/run_all_tests.sh

# 3. 查看汇总结果
bash scripts/evaluate_all_results.sh
cat results/summary_report.txt
less results/detailed_report.txt

# 4. 查看资源使用情况
ls -lh logs/
cat logs/*_resource_usage.txt | grep "Maximum resident set size"
```

### 主要特点

- 🧬 **多组学整合**：结合突变、表达和网络数据
- 🔗 **网络驱动**：利用 PPI、GRN 和互斥网络
- 🎯 **双层预测**：支持队列级和个体级分析
- 🌐 **泛癌能力**：跨多种癌症类型识别驱动基因

### 联系方式

如有问题，请通过 [GitHub Issues](https://github.com/xqxu921/HyperNetWalk/issues) 联系我们。

---

<div align="center">

**⭐ 如果 HyperNetWalk 对您有帮助，请给我们一个 Star！⭐**

Made with ❤️ by the HyperNetWalk Team

</div>