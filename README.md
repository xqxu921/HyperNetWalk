# 数据下载说明

## 自动下载（推荐）

运行下载脚本：
```bash
bash download_data.sh
```

## 网络问题备选方案

如果遇到网络问题（无法访问 AWS S3 或其他数据源），请直接下载预处理数据：

**SourceForge 镜像：**  
https://sourceforge.net/projects/hypernetwork/files/data/

下载 `data.tar.gz` 后解压到项目根目录：
```bash
tar -xzf data.tar.gz
```

## 预期目录结构

data/
├── metadata/
├── raw/
├── NETWORK/
└── DRIVER/

# 环境安装

## 克隆项目
git clone https://github.com/xqxu921/HyperNetWalk.git
cd HyperNetWalk

## 安装环境
项目包含两套环境：

Python 环境 (WeSME)

R 环境 (HyperNetWalk)，并使用 renv 管理 R 包

运行以下命令一键创建：
bash setup_environment.sh

脚本会按顺序：

创建 Python 环境 wesme

创建 R 环境 hypernetwalk

使用 renv 恢复 R 包依赖

## 激活环境
conda activate hypernetwalk

## 注意事项

需要联网以安装 Conda 包和 CRAN/Bioconductor 包。

建议 Conda 使用 conda-forge 作为主渠道以保证 R 4.5.2 的兼容性。

如果中途报错，可尝试删除已有环境再重建：
conda remove -n hypernetwalk --all
conda remove -n wesme --all
bash setup_environment.sh
