# 🔬 HyperNetWalk

**HyperNetWalk** is a hypergraph-based framework for **pan-cancer driver gene identification** across multi-omics layers.  
It integrates **protein–protein interaction (PPI)**, **gene regulatory (GRN)**, and **mutual exclusivity (ME)** networks to identify driver genes at both **cohort** and **individual** levels.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R >= 4.5.2](https://img.shields.io/badge/R-%3E%3D4.5.2-green.svg)](https://cran.r-project.org/)
[![Conda](https://img.shields.io/badge/environment-conda-orange.svg)](https://docs.conda.io/)
[![WeSME](https://img.shields.io/badge/Python-WeSME-yellow.svg)](https://sourceforge.net/projects/wesme/)

---

## 📋 Table of Contents

- [1. Environment Setup](#1-environment-setup)
- [2. Data Preparation](#2-data-preparation)
- [3. Running HyperNetWalk](#3-running-hypernetwalk)
- [4. Evaluation](#4-evaluation)
- [5. Quick Validation (Recommended)](#5-quick-validation-recommended)
- [6. Repository Structure](#6-repository-structure)
- [7. Citation](#7-citation)
- [8. 中文说明（简要）](#8-中文说明简要)

---

## 🧩 1. Environment Setup

Clone this repository and automatically configure all environments:

```bash
git clone https://github.com/xqxu921/HyperNetWalk.git
cd HyperNetWalk
bash scripts/setup_environment.sh
