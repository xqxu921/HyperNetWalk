#!/bin/bash

################################################################################
# HyperNetWalk Results Evaluation Pipeline
# 该脚本评估所有预测结果并生成汇总报告
################################################################################

set -e
set -u

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# 配置参数
RESULTS_DIR="results"
BENCHMARK="CGC"  # 可选: CGC, IntOGen, Combined
OUTPUT_SUMMARY="$RESULTS_DIR/summary_report.txt"
OUTPUT_DETAILED="$RESULTS_DIR/detailed_report.txt"

# 定义12个癌症类型
CANCER_TYPES=("BRCA" "COAD" "HNSC" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "PRAD" "STAD" "THCA" "UCEC")

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
    print_header "检查评估环境"
    
    # 检查conda环境
    if [[ "$CONDA_DEFAULT_ENV" != "hypernetwalk" ]]; then
        print_error "未激活hypernetwalk环境"
        echo "请运行: conda activate hypernetwalk"
        exit 1
    fi
    print_success "Conda环境正确: $CONDA_DEFAULT_ENV"
    
    # 检查结果目录
    if [ ! -d "$RESULTS_DIR" ]; then
        print_error "结果目录不存在: $RESULTS_DIR"
        echo "请先运行 run_all_tests.sh"
        exit 1
    fi
    print_success "结果目录存在"
    
    # 检查评估脚本
    if [ ! -f "src/evaluation.R" ]; then
        print_error "评估脚本不存在: src/evaluation.R"
        exit 1
    fi
    print_success "评估脚本检查通过"
}

evaluate_pancancer() {
    print_header "评估Pan-cancer群体结果"
    
    if [ ! -d "$RESULTS_DIR/PANCAN" ]; then
        print_error "Pan-cancer结果不存在"
        return 1
    fi
    
    Rscript src/evaluation.R \
        --mode pancancer \
        --level cohort \
        --predicted "$RESULTS_DIR/PANCAN" \
        --benchmark "$BENCHMARK" \
        --output "$RESULTS_DIR/PANCAN/evaluation_results.txt"
    
    if [ $? -eq 0 ]; then
        print_success "Pan-cancer评估完成"
        return 0
    else
        print_error "Pan-cancer评估失败"
        return 1
    fi
}

evaluate_cancer_cohort() {
    local cancer_type=$1
    print_info "评估 ${cancer_type} 群体结果..."
    
    if [ ! -d "$RESULTS_DIR/$cancer_type" ]; then
        print_error "${cancer_type} 结果目录不存在"
        return 1
    fi
    
    Rscript src/evaluation.R \
        --mode single_cancer \
        --level cohort \
        --cancer_type "$cancer_type" \
        --predicted "$RESULTS_DIR/$cancer_type/" \
        --benchmark "$BENCHMARK" \
        --output "$RESULTS_DIR/$cancer_type/cohort_evaluation_results.txt" 2>&1 | grep -v "^>"
    
    if [ $? -eq 0 ]; then
        print_success "${cancer_type} 群体评估完成"
        return 0
    else
        print_error "${cancer_type} 群体评估失败"
        return 1
    fi
}

evaluate_cancer_individual() {
    local cancer_type=$1
    print_info "评估 ${cancer_type} 个体结果..."
    
    if [ ! -d "$RESULTS_DIR/$cancer_type" ]; then
        print_error "${cancer_type} 结果目录不存在"
        return 1
    fi
    
    Rscript src/evaluation.R \
        --mode single_cancer \
        --level individual \
        --cancer_type "$cancer_type" \
        --predicted "$RESULTS_DIR/$cancer_type/" \
        --benchmark "$BENCHMARK" \
        --output "$RESULTS_DIR/$cancer_type/individual_evaluation_results.txt" 2>&1 | grep -v "^>"
    
    if [ $? -eq 0 ]; then
        print_success "${cancer_type} 个体评估完成"
        return 0
    else
        print_error "${cancer_type} 个体评估失败"
        return 1
    fi
}

extract_metrics() {
    local eval_file=$1
    local cancer=$2
    local level=$3
    
    if [ ! -f "$eval_file" ]; then
        echo "N/A"
        return
    fi
    
    # 提取关键指标（根据实际evaluation.R输出格式调整）
    # 这里是示例，需要根据实际输出格式修改
    local metrics=$(grep -E "AUC|Precision|Recall|F1" "$eval_file" | head -n 4 | awk '{print $NF}' | paste -sd ',' -)
    
    if [ -z "$metrics" ]; then
        echo "N/A"
    else
        echo "$metrics"
    fi
}

generate_summary_report() {
    print_header "生成汇总报告"
    
    {
        echo "╔════════════════════════════════════════════════════════════════╗"
        echo "║          HyperNetWalk 结果评估汇总报告                        ║"
        echo "╚════════════════════════════════════════════════════════════════╝"
        echo ""
        echo "生成时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "基准数据集: $BENCHMARK"
        echo ""
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "1. Pan-cancer 群体水平评估"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        
        if [ -f "$RESULTS_DIR/PANCAN/evaluation_results.txt" ]; then
            cat "$RESULTS_DIR/PANCAN/evaluation_results.txt" | head -n 20
        else
            echo "   [评估结果不存在]"
        fi
        
        echo ""
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "2. 各癌症类型评估结果"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""
        printf "%-10s %-20s %-20s\n" "癌症类型" "群体评估状态" "个体评估状态"
        printf "%-10s %-20s %-20s\n" "--------" "----------------" "----------------"
        
        for cancer in "${CANCER_TYPES[@]}"; do
            local cohort_status="✗ 未找到"
            local individual_status="✗ 未找到"
            
            if [ -f "$RESULTS_DIR/$cancer/cohort_evaluation_results.txt" ]; then
                cohort_status="✓ 完成"
            fi
            
            if [ -f "$RESULTS_DIR/$cancer/individual_evaluation_results.txt" ]; then
                individual_status="✓ 完成"
            fi
            
            printf "%-10s %-20s %-20s\n" "$cancer" "$cohort_status" "$individual_status"
        done
        
        echo ""
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "3. 统计摘要"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        
        local total_cancers=${#CANCER_TYPES[@]}
        local cohort_completed=0
        local individual_completed=0
        
        for cancer in "${CANCER_TYPES[@]}"; do
            [ -f "$RESULTS_DIR/$cancer/cohort_evaluation_results.txt" ] && ((cohort_completed++)) || true
            [ -f "$RESULTS_DIR/$cancer/individual_evaluation_results.txt" ] && ((individual_completed++)) || true
        done
        
        echo ""
        echo "总癌症类型数: $total_cancers"
        echo "完成群体评估: $cohort_completed / $total_cancers"
        echo "完成个体评估: $individual_completed / $total_cancers"
        echo ""
        
        if [ -d "logs" ]; then
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "4. 计算资源使用情况"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo ""
            
            if [ -f "logs/run_summary.txt" ]; then
                cat "logs/run_summary.txt" | tail -n 20
            else
                echo "   [运行日志不存在]"
            fi
        fi
        
        echo ""
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "5. 文件位置"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""
        echo "• 汇总报告: $OUTPUT_SUMMARY"
        echo "• 详细报告: $OUTPUT_DETAILED"
        echo "• 各癌症类型详细结果: $RESULTS_DIR/<CANCER_TYPE>/evaluation_results.txt"
        echo "• 运行日志: logs/"
        echo ""
        echo "╔════════════════════════════════════════════════════════════════╗"
        echo "║                     报告生成完毕                               ║"
        echo "╚════════════════════════════════════════════════════════════════╝"
        
    } > "$OUTPUT_SUMMARY"
    
    print_success "汇总报告已生成: $OUTPUT_SUMMARY"
}

generate_detailed_report() {
    print_header "生成详细报告"
    
    {
        echo "HyperNetWalk 详细评估报告"
        echo "========================================"
        echo "生成时间: $(date)"
        echo ""
        
        echo "========================================"
        echo "Pan-cancer 群体评估"
        echo "========================================"
        if [ -f "$RESULTS_DIR/PANCAN/evaluation_results.txt" ]; then
            cat "$RESULTS_DIR/PANCAN/evaluation_results.txt"
        else
            echo "评估结果不存在"
        fi
        echo ""
        
        for cancer in "${CANCER_TYPES[@]}"; do
            echo "========================================"
            echo "$cancer - 群体评估"
            echo "========================================"
            if [ -f "$RESULTS_DIR/$cancer/cohort_evaluation_results.txt" ]; then
                cat "$RESULTS_DIR/$cancer/cohort_evaluation_results.txt"
            else
                echo "评估结果不存在"
            fi
            echo ""
            
            echo "========================================"
            echo "$cancer - 个体评估"
            echo "========================================"
            if [ -f "$RESULTS_DIR/$cancer/individual_evaluation_results.txt" ]; then
                cat "$RESULTS_DIR/$cancer/individual_evaluation_results.txt"
            else
                echo "评估结果不存在"
            fi
            echo ""
        done
        
    } > "$OUTPUT_DETAILED"
    
    print_success "详细报告已生成: $OUTPUT_DETAILED"
}

compare_with_baselines() {
    print_header "与基线方法比较（可选）"
    
    # 如果有基线方法的结果，在这里添加比较逻辑
    print_info "如需与其他方法比较，请手动添加比较脚本"
}

################################################################################
# 主流程
################################################################################

main() {
    local start_time=$(date +%s)
    
    print_header "HyperNetWalk 结果评估流程"
    echo "开始时间: $(date)"
    echo ""
    
    # 检查环境
    check_prerequisites
    
    # 评估Pan-cancer
    evaluate_pancancer || print_error "Pan-cancer评估失败，继续其他评估..."
    
    # 评估各癌症类型 - 群体水平
    print_header "评估癌症类型群体结果 (${#CANCER_TYPES[@]}个)"
    local cohort_success=0
    for cancer in "${CANCER_TYPES[@]}"; do
        evaluate_cancer_cohort "$cancer"
    done
    print_info "群体评估完成"
    
    # 评估各癌症类型 - 个体水平
    print_header "评估癌症类型个体结果 (${#CANCER_TYPES[@]}个)"
    local individual_success=0
    for cancer in "${CANCER_TYPES[@]}"; do
        evaluate_cancer_individual "$cancer"
    done
    print_info "个体评估完成"
    
    # 生成报告
    generate_summary_report
    generate_detailed_report
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    print_header "评估完成！"
    print_success "总用时: ${duration}秒"
    echo ""
    print_info "查看汇总报告: cat $OUTPUT_SUMMARY"
    print_info "查看详细报告: less $OUTPUT_DETAILED"
}

# 捕获中断信号
trap 'echo -e "\n${RED}评估被中断${NC}"; exit 130' INT TERM

# 运行主程序
main "$@"