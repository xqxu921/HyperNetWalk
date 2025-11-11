#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 数据下载脚本 - TCGA/PANCAN 癌症基因组学数据
# ============================================================
# 
# 【网络问题备选方案】
# 如果遇到网络连接问题（如无法访问 AWS S3 或 StringDB），
# 可以直接从 SourceForge 下载已处理的数据：
# 
#   https://sourceforge.net/projects/hypernetwork/files/data/processed_data/
# 
# 下载后解压到项目根目录即可，无需运行此脚本。
# ============================================================

# ---------------- 配置 ----------------
BASE_DIR="data"
RAW_DIR="$BASE_DIR/raw"
NETWORK_DIR="$BASE_DIR/NETWORK"
DRIVER_DIR="$BASE_DIR/DRIVER"
GENE_MAP_DIR="$BASE_DIR/metadata"

MAX_TRIES=6
INITIAL_SLEEP=3
ARIA2C_SEGMENTS=8

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[1;34m'; NC='\033[0m'

mkdir -p "$RAW_DIR" "$NETWORK_DIR" "$DRIVER_DIR" "$GENE_MAP_DIR" "$BASE_DIR/logs"
ERROR_LOG="$BASE_DIR/logs/download_errors.log"
: > "$ERROR_LOG"

has_cmd() { command -v "$1" >/dev/null 2>&1; }

echo -e "${BLUE}检测下载工具...${NC}"
if has_cmd aria2c; then
    DL_TOOL="aria2c"; echo -e "${GREEN}使用 aria2c${NC}"
elif has_cmd curl; then
    DL_TOOL="curl"; echo -e "${GREEN}使用 curl${NC}"
elif has_cmd wget; then
    DL_TOOL="wget"; echo -e "${GREEN}使用 wget${NC}"
else
    echo -e "${YELLOW}未检测到下载工具，建议安装 aria2c/curl/wget${NC}"; exit 1
fi

# ---------------- 下载函数 ----------------
download_one() {
    local url="$1"
    local out="$2"
    local desc="$3"
    mkdir -p "$(dirname "$out")"
    if [ -f "$out" ] && [ -s "$out" ]; then
        echo -e "${BLUE}[跳过]${NC} 已存在: $out"; return 0
    fi

    local try=1
    local sleep_time=$INITIAL_SLEEP

    while [ $try -le $MAX_TRIES ]; do
        echo -e "${YELLOW}[$try/$MAX_TRIES] 下载: $desc${NC}"
        if has_cmd aria2c; then
            aria2c --file-allocation=none -x"$ARIA2C_SEGMENTS" -s"$ARIA2C_SEGMENTS" -d "$(dirname "$out")" -o "$(basename "$out")" --max-tries=1 --retry-wait=5 "$url" >/dev/null 2>&1 || true
        fi
        if [ ! -s "$out" ] && has_cmd curl; then
            curl --fail --location --connect-timeout 60 --max-time 600 -o "$out" "$url" || true
        fi
        if [ ! -s "$out" ] && has_cmd wget; then
            wget -c -O "$out" "$url" || true
        fi
        if [ -s "$out" ]; then
            echo -e "${GREEN}✓ 下载成功: $desc -> $out${NC}"; return 0
        fi
        # IPv4 回退
        echo -e "${YELLOW}尝试 IPv4 回退...${NC}"
        if has_cmd curl; then curl --ipv4 --fail -o "$out" "$url" || true; fi
        if [ ! -s "$out" ] && has_cmd wget; then wget -4 -c -O "$out" "$url" || true; fi
        if [ -s "$out" ]; then echo -e "${GREEN}✓ (IPv4) 下载成功: $desc -> $out${NC}"; return 0; fi
        echo -e "${YELLOW}下载失败（尝试 $try），等待 ${sleep_time}s 后重试...${NC}"
        sleep $sleep_time; sleep_time=$((sleep_time*2)); try=$((try+1))
    done

    echo -e "${RED}✗ 最终失败: $desc${NC}"
    echo -e "FAILED\t$desc\t$url" >> "$ERROR_LOG"
    return 1
}

# ---------------- 解压函数 ----------------
decompress_gz() {
    local file="$1"
    local desc="$2"
    if [ -f "$file" ] && [ -s "$file" ]; then
        echo -e "${YELLOW}解压: $desc${NC}"
        gunzip -kf "$file" 2>/dev/null || true
        # 检查解压是否成功
        local uncompressed="${file%.gz}"
        if [ -f "$uncompressed" ] && [ -s "$uncompressed" ]; then
            echo -e "${GREEN}✓ 解压成功: $uncompressed${NC}"
        else
            echo -e "${RED}✗ 解压失败: $file${NC}"
            echo -e "DECOMPRESS_FAILED\t$desc\t$file" >> "$ERROR_LOG"
        fi
    fi
}

# ---------------- 下载 ----------------
echo -e "${BLUE}开始批量下载...${NC}"

# gene mapping
download_one "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap" \
    "$GENE_MAP_DIR/gencode.v36.annotation.gtf.gene.probemap" "GENCODE probeMap" || true

# TCGA & PANCAN
CANCERS=(BRCA COAD HNSC KIRC KIRP LIHC LUAD LUSC PRAD STAD THCA UCEC)
for c in "${CANCERS[@]}"; do
    download_one "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-${c}.star_counts.tsv.gz" \
        "$RAW_DIR/${c}.star_counts.tsv.gz" "$c counts" || true
    decompress_gz "$RAW_DIR/${c}.star_counts.tsv.gz" "$c counts"
    
    download_one "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-${c}.star_tpm.tsv.gz" \
        "$RAW_DIR/${c}.star_tpm.tsv.gz" "$c TPM" || true
    decompress_gz "$RAW_DIR/${c}.star_tpm.tsv.gz" "$c TPM"
    
    download_one "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-${c}.somaticmutation_wxs.tsv.gz" \
        "$RAW_DIR/${c}.somaticmutation_wxs.tsv.gz" "$c mutation" || true
    decompress_gz "$RAW_DIR/${c}.somaticmutation_wxs.tsv.gz" "$c mutation"
    
    download_one "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-${c}.clinical.tsv.gz" \
        "$RAW_DIR/${c}.clinical.tsv.gz" "$c clinical" || true
    decompress_gz "$RAW_DIR/${c}.clinical.tsv.gz" "$c clinical"
done

# PANCAN
download_one "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz" \
    "$RAW_DIR/EB_AdjustPANCAN_RNASeqV2.geneExp.xena.gz" "PANCAN expression" || true
decompress_gz "$RAW_DIR/EB_AdjustPANCAN_RNASeqV2.geneExp.xena.gz" "PANCAN expression"

download_one "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz" \
    "$RAW_DIR/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz" "PANCAN mutation" || true
decompress_gz "$RAW_DIR/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz" "PANCAN mutation"

# NETWORK
download_one "https://stringdb-downloads.org/download/stream/protein.links.v12.0/9606.protein.links.v12.0.onlyAB.tsv.gz" \
    "$NETWORK_DIR/9606.protein.links.v12.0.onlyAB.tsv.gz" "STRING PPI" || true
decompress_gz "$NETWORK_DIR/9606.protein.links.v12.0.onlyAB.tsv.gz" "STRING PPI"

download_one "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz" \
    "$GENE_MAP_DIR/9606.protein.info.v12.0.txt.gz" "STRING protein info" || true
decompress_gz "$GENE_MAP_DIR/9606.protein.info.v12.0.txt.gz" "STRING protein info"

download_one "http://www.zpliulab.cn/RegNetwork/static/data/human_TF_Target.7z" \
    "$NETWORK_DIR/human_TF_Target.7z" "RegNetwork" || true
if [ -f "$NETWORK_DIR/human_TF_Target.7z" ] && has_cmd 7z; then
    echo -e "${YELLOW}解压 RegNetwork...${NC}"
    7z x -y "$NETWORK_DIR/human_TF_Target.7z" -o"$NETWORK_DIR" >/dev/null 2>&1
    if [ -f "$NETWORK_DIR/human_TF_Target.txt" ]; then
        echo -e "${GREEN}✓ RegNetwork 解压成功${NC}"
    else
        echo -e "${RED}✗ RegNetwork 解压失败${NC}"
    fi
fi

# DRIVER
download_one "https://cancer.sanger.ac.uk/cosmic/census/all?export=tsv" \
    "$DRIVER_DIR/CGC_Tier1.tsv" "COSMIC CGC Tier1" || true

download_one "https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip" \
    "$DRIVER_DIR/IntOGen-Drivers-20240920.zip" "IntOGen Drivers" || true

if [ -f "$DRIVER_DIR/IntOGen-Drivers-20240920.zip" ] && has_cmd unzip; then
    echo -e "${YELLOW}解压 IntOGen...${NC}"
    unzip -j -o "$DRIVER_DIR/IntOGen-Drivers-20240920.zip" \
        "2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv" \
        -d "$DRIVER_DIR" >/dev/null 2>&1
    if [ -f "$DRIVER_DIR/Compendium_Cancer_Genes.tsv" ]; then
        echo -e "${GREEN}✓ IntOGen 解压成功${NC}"
    else
        echo -e "${RED}✗ IntOGen 解压失败${NC}"
    fi
fi

# ---------------- 完成 ----------------
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}下载流程结束${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# 检查是否有失败项
if [ -s "$ERROR_LOG" ]; then
    echo -e "${RED}⚠ 存在失败项，请检查: $ERROR_LOG${NC}"
    echo -e "${YELLOW}如果网络问题无法解决，可从以下地址下载完整数据：${NC}"
    echo -e "${BLUE}https://sourceforge.net/projects/hypernetwork/files/data/processed_data/${NC}"
else
    echo -e "${GREEN}✓ 所有文件下载成功！${NC}"
fi

echo ""
echo -e "${BLUE}最终目录结构检查：${NC}"
tree -L 2 "$BASE_DIR" 2>/dev/null || find "$BASE_DIR" -type f | sort