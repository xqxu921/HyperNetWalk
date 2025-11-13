#!/bin/bash

################################################################################
# HyperNetWalk Complete Test Pipeline
# 该脚本依次运行Pan-cancer和12个癌症类型的群体及个体预测
################################################################################

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 配置参数
CORES=64
INPUT_DIR="data/processed"
PPI_NETWORK="data/NETWORK/STRINGv12.txt"
GRN_NETWORK="data/NETWORK/RegNet_human_V2.txt"
OUTPUT_DIR="results"
LOG_DIR="logs/STRING_v12.txt"
GRN_NETWORK="data/NETWORK/RegNet_human_V2.txt"
OUTPUT_DIR="results"
LOG_DIR="logs"

# 定义12个癌症类型
CANCER_TYPES=("BRCA" "COAD" "HNSC" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "PRAD" "STAD" "THCA" "UCEC")

# 检测系统中的time命令
if [ -f /usr/bin/time ]; then
    TIME_CMD="/usr/bin/time -v"
elif command -v gtime &> /dev/null; then
    TIME_CMD="gtime -v"  # macOS上通过brew安装的GNU time
else
    echo -e "${YELLOW}Warning: /usr/bin/time not found, using bash builtin time (no memory stats)${NC}"
    TIME_CMD="time"
fi

################################################################################
# 函数定义
################################################################################

print_header() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo ""
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_info() {
    echo -e "${YELLOW}ℹ $1${NC}"
}

check_prerequisites() {
    print_header "检查运行环境"
    
    # 检查conda环境
    if [[ "$CONDA_DEFAULT_ENV" != "hypernetwalk" ]]; then
        print_error "未激活hypernetwalk环境"
        echo "请运行: conda activate hypernetwalk"
        exit 1
    fi
    print_success "Conda环境正确: $CONDA_DEFAULT_ENV"
    
    # 检查R
    if ! command -v R &> /dev/null; then
        print_error "未找到R"
        exit 1
    fi
    print_success "R版本: $(R --version | head -n1)"
    
    # 检查必要目录
    for dir in "$INPUT_DIR" "data/NETWORK" "src"; do
        if [ ! -d "$dir" ]; then
            print_error "目录不存在: $dir"
            exit 1
        fi
    done
    print_success "必要目录检查通过"
    
    # 检查网络文件
    if [ ! -f "$PPI_NETWORK" ]; then
        print_error "PPI网络文件不存在: $PPI_NETWORK"
        exit 1
    fi
    if [ ! -f "$GRN_NETWORK" ]; then
        print_error "GRN网络文件不存在: $GRN_NETWORK"
        exit 1
    fi
    print_success "网络文件检查通过"
    
    # 创建输出目录
    mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
    print_success "输出目录已创建"
}

run_pancancer_cohort() {
    print_header "运行Pan-cancer群体水平预测"
    
    local start_time=$(date +%s)
    local log_file="$LOG_DIR/pancan_cohort_resource.txt"
    
    $TIME_CMD -o "$log_file" \
    Rscript src/run_hypernetwalk.R \
        --mode pancancer \
        --level cohort \
        --input "$INPUT_DIR" \
        --ppi "$PPI_NETWORK" \
        --grn "$GRN_NETWORK" \
        --output "$OUTPUT_DIR/" \
        --cores "$CORES" 2>&1 | tee "$LOG_DIR/pancan_cohort.log"
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    print_success "Pan-cancer群体预测完成 (用时: ${duration}秒)"
    
    # 提取内存使用信息
    if [ -f "$log_file" ] && grep -q "Maximum resident set size" "$log_file"; then
        local max_mem=$(grep "Maximum resident set size" "$log_file" | awk '{print $6}')
        print_info "峰值内存使用: $((max_mem / 1024)) MB"
    fi
}
run_cancer_type_all() {
    local cancer_type=$1
    print_header "运行 ${cancer_type} 群体水平预测"
    
    local start_time=$(date +%s)
    local log_file="$LOG_DIR/${cancer_type}_cohort_resource.txt"
    
    $TIME_CMD -o "$log_file" \
    Rscript src/run_hypernetwalk.R \
        --mode single_cancer \
        --level all \
        --cancer_type "$cancer_type" \
        --input "$INPUT_DIR/" \
        --ppi "$PPI_NETWORK" \
        --grn "$GRN_NETWORK" \
        --output "$OUTPUT_DIR/" \
        --cores "$CORES" 2>&1 | tee "$LOG_DIR/${cancer_type}_cohort&individual.log"
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    print_success "${cancer_type} 群体+个体预测完成 (用时: ${duration}秒)"
    
    if [ -f "$log_file" ] && grep -q "Maximum resident set size" "$log_file"; then
        local max_mem=$(grep "Maximum resident set size" "$log_file" | awk '{print $6}')
        print_info "峰值内存使用: $((max_mem / 1024)) MB"
    fi
}

run_cancer_type_cohort() {
    local cancer_type=$1
    print_header "运行 ${cancer_type} 群体水平预测"
    
    local start_time=$(date +%s)
    local log_file="$LOG_DIR/${cancer_type}_cohort_resource.txt"
    
    $TIME_CMD -o "$log_file" \
    Rscript src/run_hypernetwalk.R \
        --mode single_cancer \
        --level cohort \
        --cancer_type "$cancer_type" \
        --input "$INPUT_DIR/" \
        --ppi "$PPI_NETWORK" \
        --grn "$GRN_NETWORK" \
        --output "$OUTPUT_DIR/" \
        --cores "$CORES" 2>&1 | tee "$LOG_DIR/${cancer_type}_cohort.log"
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    print_success "${cancer_type} 群体预测完成 (用时: ${duration}秒)"
    
    if [ -f "$log_file" ] && grep -q "Maximum resident set size" "$log_file"; then
        local max_mem=$(grep "Maximum resident set size" "$log_file" | awk '{print $6}')
        print_info "峰值内存使用: $((max_mem / 1024)) MB"
    fi
}

run_cancer_type_individual() {
    local cancer_type=$1
    print_header "运行 ${cancer_type} 个体水平预测"
    
    local start_time=$(date +%s)
    local log_file="$LOG_DIR/${cancer_type}_individual_resource.txt"
    
    $TIME_CMD -o "$log_file" \
    Rscript src/run_hypernetwalk.R \
        --mode single_cancer \
        --level individual \
        --cancer_type "$cancer_type" \
        --input "$INPUT_DIR/" \
        --ppi "$PPI_NETWORK" \
        --grn "$GRN_NETWORK" \
        --output "$OUTPUT_DIR/" \
        --cores "$CORES" 2>&1 | tee "$LOG_DIR/${cancer_type}_individual.log"
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    print_success "${cancer_type} 个体预测完成 (用时: ${duration}秒)"
    
    if [ -f "$log_file" ] && grep -q "Maximum resident set size" "$log_file"; then
        local max_mem=$(grep "Maximum resident set size" "$log_file" | awk '{print $6}')
        print_info "峰值内存使用: $((max_mem / 1024)) MB"
    fi
}

generate_summary() {
    print_header "生成运行摘要"
    
    local summary_file="$LOG_DIR/run_summary.txt"
    
    {
        echo "HyperNetWalk 测试运行摘要"
        echo "======================================"
        echo "运行时间: $(date)"
        echo "使用核心数: $CORES"
        echo ""
        echo "完成的任务:"
        echo "  - Pan-cancer 群体预测: ✓"
        echo "  - 癌症类型数量: ${#CANCER_TYPES[@]}"
        echo "  - 群体预测: ${#CANCER_TYPES[@]} 个癌症类型"
        echo "  - 个体预测: ${#CANCER_TYPES[@]} 个癌症类型"
        echo ""
        echo "资源使用统计:"
        echo "--------------------------------------"
        
        if [ -f "$LOG_DIR/pancan_cohort_resource.txt" ]; then
            echo "Pan-cancer:"
            grep "Elapsed (wall clock) time\|Maximum resident set size" "$LOG_DIR/pancan_cohort_resource.txt" || echo "  未找到资源信息"
            echo ""
        fi
        
        for cancer in "${CANCER_TYPES[@]}"; do
            echo "${cancer}:"
            if [ -f "$LOG_DIR/${cancer}_cohort_resource.txt" ]; then
                echo "  Cohort:"
                grep "Elapsed (wall clock) time\|Maximum resident set size" "$LOG_DIR/${cancer}_cohort_resource.txt" || echo "    未找到资源信息"
            fi
            if [ -f "$LOG_DIR/${cancer}_individual_resource.txt" ]; then
                echo "  Individual:"
                grep "Elapsed (wall clock) time\|Maximum resident set size" "$LOG_DIR/${cancer}_individual_resource.txt" || echo "    未找到资源信息"
            fi
            echo ""
        done
        
    } > "$summary_file"
    
    print_success "摘要已保存到: $summary_file"
    cat "$summary_file"
}

################################################################################
# 主流程
################################################################################

main() {
    local overall_start=$(date +%s)
    
    print_header "HyperNetWalk 完整测试流程"
    echo "开始时间: $(date)"
    echo "使用核心数: $CORES"
    echo "癌症类型: ${CANCER_TYPES[*]}"
    echo ""
    
    # 检查环境
    check_prerequisites
    
    # 1. Pan-cancer群体预测
    run_pancancer_cohort

    # 2. 各癌症类型群体+个体预测
    print_header "开始癌症类型群体+个体预测 (${#CANCER_TYPES[@]}个)"
    for cancer in "${CANCER_TYPES[@]}"; do
        run_cancer_type_all "$cancer"
    done
    
    # # 2. 各癌症类型群体预测
    # print_header "开始癌症类型群体预测 (${#CANCER_TYPES[@]}个)"
    # for cancer in "${CANCER_TYPES[@]}"; do
    #     run_cancer_type_cohort "$cancer"
    # done
    
    # # 3. 各癌症类型个体预测
    # print_header "开始癌症类型个体预测 (${#CANCER_TYPES[@]}个)"
    # for cancer in "${CANCER_TYPES[@]}"; do
    #     run_cancer_type_individual "$cancer"
    # done
    
    # 生成摘要
    generate_summary
    
    local overall_end=$(date +%s)
    local total_duration=$((overall_end - overall_start))
    local hours=$((total_duration / 3600))
    local minutes=$(((total_duration % 3600) / 60))
    local seconds=$((total_duration % 60))
    
    print_header "所有测试完成！"
    print_success "总用时: ${hours}小时 ${minutes}分钟 ${seconds}秒"
    print_info "请运行 'bash evaluate_all_results.sh' 来评估结果"
}

# 捕获中断信号
trap 'echo -e "\n${RED}测试被中断${NC}"; exit 130' INT TERM

# 运行主程序
main "$@"