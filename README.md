# TCGA-BRCA 临床—组学联合生存分析（练手项目）

作者：蔡璋睿（宁波大学 经济学 2024 级）｜2026.11｜数据：TCGA（经 GDC 开放层，无需认证）

入门练手项目：用 GDC API 下载 TCGA-BRCA（乳腺浸润癌）开放层数据，完成从原始表达计数到临床—组学联合生存分析的完整流程。灵感来自肿瘤方向博士师姐的建议——TCGA 是临床数据入门最好的免费数据源。

## 研究问题

在 TCGA-BRCA 原发肿瘤样本中，乳腺癌关键基因（**ESR1**、**ERBB2/HER2**、**MKI67**）的表达水平是否与患者总生存（OS）相关？在调整年龄和肿瘤分期后，关联是否仍然成立？

## 数据与方法选择

**数据**
- 临床：GDC `cases` 接口，TCGA-BRCA 全部 1,098 例（`vital_status`、`days_to_death`、`diagnoses.days_to_last_follow_up`、`follow_ups.days_to_follow_up`、`age_at_diagnosis`、`ajcc_pathologic_stage`）
- 表达：GDC `files` + `/data` 接口，100 例原发肿瘤（Primary Tumor）RNA-Seq **STAR - Counts** 文件（开放层，无需认证；本机网络受限实际成功下载 39 个）

**为什么这样选方法**
- 结果变量是"时间 + 是否发生事件（死亡）"→ 经典生存分析框架：Kaplan-Meier 估计、对数秩检验、Cox 比例风险模型，而不是普通回归（普通回归会丢弃删失信息、产生偏误）
- 表达计数存在文库大小差异 → CPM 标准化后取 log2（练手级处理；正式研究应考虑 TMM/DESeq2 标准化与批次校正）
- Cox 模型采用 **Efron 结近似**（并列死亡时间更准），Newton-Raphson 求解；实现已经过模拟数据验证（估计值回收真实参数，Wald z² 与对数秩 χ² 一致）
- 目标基因中位数分组 → 非参数 KM + 对数秩（不假设函数形式）；多变量 Cox 调整混杂（年龄、分期）——对应"混杂因素、效应量与置信区间"的岗位要求

## 本人贡献

独立完成：GDC API 数据获取（含过滤语法调试）、STAR - Counts 解析与 CPM 标准化、临床表构建（风险集口径统一）、KM/对数秩检验实现、Cox PH 模型实现（Efron 结）、结果图表与报告。

## 结果与局限（诚实披露）

实际入组：临床—表达成功匹配 **38 例**，死亡事件仅 **5 个**（晚期 8 例、其中 2 例死亡）。结果见 `results/summary.json` 与 `figures/km_genes.png`。

**主要结果（练手演示，不构成生物学结论）**
- 对数秩检验：ESR1 p=0.648，ERBB2 p=0.515，MKI67 p=0.609 —— 三个基因中位数分组均无显著生存差异（与预期一致）
- 单变量 Cox：ESR1 高表达 HR=0.53 [0.06, 5.13]，ERBB2 HR=1.37 [0.22, 8.65]，MKI67 HR=1.31 [0.22, 7.89] —— 点估计方向与文献大体一致（ESR1 高表达偏保护），但置信区间极宽、全部不显著
- 未报告多变量 Cox：事件数=5，按 EPV≥10 经验法则最多容纳 1 个协变量；强行放入年龄+分期会出现完全分离（某亚组零事件导致 HR 发散），故仅做单变量模型并在输出中说明原因——这是"知道什么时候不该报数"的一部分

局限：
1. **样本量小（38 例）且事件数极少（5 个）**，统计功效接近零，置信区间宽到没有信息量；本项目的价值在流程跑通与方法验证，不在结论
2. 仅取单癌种原发肿瘤，无训练/测试划分、无外部验证
3. CPM 标准化较粗糙，未处理批次效应；中位数分组丢失表达量信息、且二分化本身损失统计功效
4. 单基因分析存在多重检验问题（3 个基因），正式研究应做 FDR 校正
5. GDC 过滤语法中 `in` 数组形式在当前部署返回 HTTP 500，本项目使用 `=` 对象形式（已写入脚本注释）；临床随访时长需经 `expand=follow_ups` 获取，字段名是 `days_to_follow_up` 而非 `days_to_last_follow_up`

## 复现

```bash
pip install -r requirements.txt
python scripts/01_download.py   # 下载（约 420MB，逐文件断点续传）
python scripts/02_analyze.py    # 分析 → results/ 与 figures/
```

## 目录

```
scripts/01_download.py   GDC API 数据下载（临床 + 表达清单 + 逐文件下载）
scripts/02_analyze.py    生存分析全流程（KM / 对数秩 / Cox PH）
data/                    原始数据（体积大，不入 git）
results/summary.json     结构化结果
figures/km_genes.png     KM 曲线
```
