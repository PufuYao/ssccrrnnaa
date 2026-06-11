# scRNA-seq 数据处理 Pipeline

> **论文**: Single-cell transcriptome atlas of spontaneous dry age-related macular degeneration in macaques  
> **测序平台**: 10X Genomics Single Cell 3' v3, Illumina HiSeq 4000  
> **物种**: Macaca fascicularis (食蟹猴)  
> **数据**: GSA CRA010812, 5 个实验 × 4 lane = 40 个 fastq 文件

---

## Pipeline 总览

```
Fastq (GSA下载)
    │
    ▼
Step 1  Cell Ranger count ─── 比对 + 定量 → 基因-细胞矩阵
    │
    ▼
Step 2  QC 过滤 ─── 去除低质量细胞和基因
    │
    ▼
Step 3  标准化与整合 ─── LogNormalize, 批次校正
    │
    ▼
Step 4  降维与聚类 ─── PCA → UMAP → 无监督聚类
    │
    ▼
Step 5  细胞类型注释 ─── 标记基因鉴定 9 种细胞类型
    │
    ▼
Step 6  差异表达分析 ─── AMD vs Control, 各细胞类型 DEG
    │
    ▼
Step 7  GO 富集分析 ─── 功能通路注释
    │
    ▼
Step 8  SASP 评分 ─── 衰老相关分泌表型
    │
    ▼
Step 9  CellPhoneDB ─── 细胞间配体-受体通讯
    │
    ▼
Step 10 可视化与验证 ─── 小鼠模型验证
```

---

## Step 0: 数据下载与准备

### 从 GSA 下载 fastq 文件

```bash
# GSA accession: CRA010812
# 5 个实验, 每个 4 个 lane, 共 40 个 fastq.gz 文件
# 文件命名: {Run}_f1.fastq.gz (Read 1) 和 {Run}_r2.fastq.gz (Read 2)

# 使用 aspera 或 wget 下载
wget -r -np -nd https://ngdc.cncb.ac.cn/gsa/browse/CRA010812
```

### 组织文件结构 (Cell Ranger 要求)

```bash
# Cell Ranger count 要求每个样本一个目录, 内含 fastq 文件
# 命名格式: {SampleID}_S1_L00{Lane}_{ReadType}_001.fastq.gz

mkdir -p data/{N1,N2,A1,A2,A3}

# 软链接并重命名 fastq 文件
# 示例: N1 样本 (CRR771605-771608)
ln -s /path/to/CRR771605_f1.fastq.gz data/N1/N1_S1_L001_R1_001.fastq.gz
ln -s /path/to/CRR771605_r2.fastq.gz data/N1/N1_S1_L001_R2_001.fastq.gz
# ... 依次处理 L002-L004
```

---

## Step 1: Cell Ranger count — 比对与定量

### 为什么做这一步？

10X Genomics 的 fastq 文件包含三个关键信息：
- **Read 1 (R1, 26 bp)**: 16 bp 细胞条形码 (barcode) + 10 bp UMI
- **Read 2 (R2, ~90 bp)**: cDNA 插入片段（转录本序列）

Cell Ranger count 完成三项核心任务：
1. **barcode 识别与纠错**：从 ~737,000 种可能的 barcode 中识别真正包裹了细胞的 barcode
2. **STAR 比对**：将 R2 比对到参考基因组
3. **UMI 计数**：基于 UMI 去重，生成每个细胞中每个基因的 UMI 计数矩阵

### 执行

```bash
# 下载食蟹猴参考基因组和注释 (Macaca fascicularis)
# 可使用 Macaca_fascicularis_6.0 (macFas6)
wget ftp://ftp.ensembl.org/pub/release-110/fasta/macaca_fascicularis/dna/Macaca_fascicularis_Macaca_fascicularis_6.0.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-110/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_6.0.110.gtf.gz

# 构建 Cell Ranger 参考索引
cellranger mkref \
    --genome=macFas6 \
    --fasta=Macaca_fascicularis_Macaca_fascicularis_6.0.dna.toplevel.fa \
    --genes=Macaca_fascicularis.Macaca_fascicularis_6.0.110.gtf \
    --nthreads=16

# 每个样本独立运行 cellranger count
# 10X v3 化学试剂: --chemistry=SC3Pv3
# 预期捕获 ~10,000 cells/样本

for sample in N1 N2 A1 A2 A3; do
    cellranger count \
        --id=${sample} \
        --transcriptome=macFas6 \
        --fastqs=data/${sample} \
        --sample=${sample} \
        --expect-cells=10000 \
        --chemistry=SC3Pv3 \
        --localcores=16 \
        --localmem=64
done
```

### 输出

```
{sample}/outs/
├── filtered_feature_bc_matrix/   ← 过滤后的基因-细胞矩阵 (核心输出)
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── raw_feature_bc_matrix/        ← 原始矩阵 (含空液滴)
├── web_summary.html              ← QC 报告
└── possorted_genome_bam.bam      ← 比对 BAM 文件
```

`filtered_feature_bc_matrix` 中 Cell Ranger 已根据 UMI 拐点自动区分"真实细胞"和"空液滴"。

---

## Step 2: 读入数据与 QC 过滤

### 为什么做这一步？

Cell Ranger 的自动过滤并不完美，需要进一步手动 QC：
- **低质量细胞**：UMI 太少说明细胞已死或膜破裂
- **双细胞 (doublets)**：一个液滴包裹两个细胞
- **高线粒体基因比例**：>20% 线粒体 reads 说明细胞正在凋亡
- **核糖体基因污染**：技术噪音

### 使用 Seurat (R) 执行

```r
library(Seurat)
library(dplyr)

# 读入所有样本
samples <- c("N1", "N2", "A1", "A2", "A3")
seus <- list()
for (s in samples) {
    counts <- Read10X(file.path(s, "outs/filtered_feature_bc_matrix"))
    seus[[s]] <- CreateSeuratObject(counts, project = s, min.cells = 3, min.features = 200)
    seus[[s]]$sample <- s
    seus[[s]]$group <- ifelse(grepl("^A", s), "AMD", "Control")
}

# 合并所有样本
merged <- merge(seus[[1]], seus[2:5])

# 计算 QC 指标
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged[["percent.ribo"]] <- PercentageFeatureSet(merged, pattern = "^RP[SL]")

# 可视化 QC
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "sample", pt.size = 0.1)

# 过滤标准 (根据质控图调整)
merged <- subset(merged, 
    nFeature_RNA > 500 &          # 最少 500 个基因
    nFeature_RNA < 4000 &         # 去除可能的 doublets
    nCount_RNA < 20000 &          # 去除高 UMI 异常值
    percent.mt < 15                # 线粒体比例 < 15%
)
```

论文最终保留约 **30,000 个高质量视网膜细胞** 进行下游分析。

---

## Step 3: 标准化与样本整合

### 为什么做这一步？

- **标准化**：不同细胞的测序深度不同（UMI 总数差异大），需要对文库大小进行归一化，使细胞间可比
- **LogNormalize**：对数转换使数据更接近正态分布，放大小差异信号（论文使用的标准化方法）
- **样本整合**：消除不同样本间的技术批次效应，保留真实的生物学差异

```r
# 拆分为样本列表进行整合
obj.list <- SplitObject(merged, split.by = "sample")

# 每个样本独立标准化
obj.list <- lapply(obj.list, function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    return(x)
})

# 选择用于整合的特征 (anchors)
features <- SelectIntegrationFeatures(object.list = obj.list)

# 寻找整合锚点并进行 CCA 校正
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(integrated) <- "integrated"
```

> **为什么使用 CCA (Canonical Correlation Analysis) 整合？**  
> CCA 识别不同数据集中共享的相关结构，比简单合并更能保留各样本共同的生物学变异，同时消除技术性批次效应。论文中提及未观察到显著的批次效应，但仍建议进行整合。

---

## Step 4: 降维与聚类

### 为什么做这一步？

- **降维**：scRNA-seq 数据有 ~20,000 个基因，是"高维诅咒"问题。PCA 降维到少量主成分保留主要变异信号
- **UMAP**：将高维空间中的细胞距离关系映射到二维平面，可视化细胞分布
- **聚类**：基于表达相似性将细胞分组，发现细胞类型

```r
# 数据缩放 (去除细胞间总表达量差异)
integrated <- ScaleData(integrated, vars.to.regress = c("nCount_RNA", "percent.mt"))

# PCA 降维
integrated <- RunPCA(integrated, npcs = 50)

# 选择显著的主成分 (基于 elbow plot 或 JackStraw)
ElbowPlot(integrated, ndims = 50)
# 通常选择 20-30 个 PC

# UMAP 可视化
integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)

# 无监督聚类 (Louvain 算法)
integrated <- FindClusters(integrated, resolution = 0.8)

# 可视化
DimPlot(integrated, group.by = "seurat_clusters", label = TRUE)
DimPlot(integrated, group.by = "group")  # AMD vs Control
```

---

## Step 5: 细胞类型注释

### 为什么做这一步？

聚类给出的是无意义的数字编号 (cluster 0, 1, 2...)，需要根据已知的视网膜细胞标记基因将每个 cluster 注释为生物学上有意义的细胞类型。

### 论文鉴定的 9 种细胞类型及其标记基因

```r
# 找出每个 cluster 的差异高表达基因
markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 视网膜细胞类型标记基因
celltype_markers <- list(
    "Rod"          = c("RHO", "GNAT1", "PDE6A", "CNGA1"),    # 视杆细胞
    "Cone"         = c("OPN1SW", "OPN1LW", "ARR3", "PDE6H"), # 视锥细胞
    "Bipolar"      = c("VSX2", "TRPM1", "GRM6", "CABP5"),    # 双极细胞
    "Amacrine"     = c("GAD1", "GAD2", "SLC6A9", "TFAP2A"),  # 无长突细胞
    "Horizontal"   = c("LHX1", "CALB1", "ONECUT1"),          # 水平细胞
    "RGC"          = c("RBPMS", "THY1", "NEFL", "POU4F2"),   # 视网膜神经节细胞
    "Muller"       = c("RLBP1", "SLC1A3", "GLUL", "APOE"),   # Müller 胶质细胞
    "Microglia"    = c("CX3CR1", "P2RY12", "TREM2", "C1QA"), # 小胶质细胞
    "RPE"          = c("RPE65", "LRAT", "BEST1", "OTX2")     # 视网膜色素上皮
)

# 注释
new_ids <- c("0" = "Rod", "1" = "Cone", ...)  # 根据标记基因匹配
integrated$celltype <- RenameIdents(integrated, new_ids)
```

> **关键点**：注释需要结合文献已知的视网膜细胞标记基因和 cluster 的差异高表达基因进行人工判断，不是纯自动化的过程。

---

## Step 6: 差异表达分析 (DEG)

### 为什么做这一步？

论文的核心目标：**比较干性 AMD 组与正常对照组在每个细胞类型中的转录组差异**。

```r
# 对每个细胞类型做 AMD vs Control 差异表达
Idents(integrated) <- "celltype"
degs_list <- list()

for (ct in unique(integrated$celltype)) {
    sub <- subset(integrated, celltype == ct)
    Idents(sub) <- "group"
    
    degs <- FindMarkers(sub, 
        ident.1 = "AMD", 
        ident.2 = "Control",
        min.pct = 0.1,
        logfc.threshold = 0.25,
        test.use = "wilcox"      # Wilcoxon 秩和检验
    )
    degs$celltype <- ct
    degs$gene <- rownames(degs)
    degs_list[[ct]] <- degs
}

all_degs <- bind_rows(degs_list)

# 筛选显著 DEG: |log2FC| > 0.25 且 adjusted p < 0.05
sig_degs <- all_degs %>% 
    filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)

# 统计各细胞类型的 DEG 数量
sig_degs %>% 
    group_by(celltype) %>% 
    summarise(up = sum(avg_log2FC > 0), down = sum(avg_log2FC < 0))
```

### 论文的主要 DEG 发现

| 细胞类型 | 上调 DEG 数 | 下调 DEG 数 | 关键基因 |
|----------|------------|------------|----------|
| BCs | 多 | 多 | — |
| Rod | 多 | — | 光敏感相关基因 |
| Müller glia | 多 | — | KDR (VEGFR2) |
| HCs | — | 多 | — |
| RGCs | — | 多 | — |

---

## Step 7: GO 富集分析

### 为什么做这一步？

DEG 是一长串基因列表，需要归纳为更高层次的功能通路，理解 AMD 背后的生物学过程。

```r
library(clusterProfiler)
library(org.Mmu.eg.db)  # 食蟹猴需转换到同源基因, 或用 human org.Hs.eg.db

# 功能: 对每个细胞类型的上/下调 DEG 分别做 GO 富集
do_go <- function(genes, ontology = "BP") {
    enrichGO(
        gene = genes,
        OrgDb = org.Hs.eg.db,  # 使用人类同源注释
        keyType = "SYMBOL",
        ont = ontology,        # BP=生物学过程, MF=分子功能, CC=细胞组分
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1
    )
}

# RPE 上调基因 GO
rpe_up_go <- do_go(rpe_up_genes)
dotplot(rpe_up_go)  # 显示 Angiogenesis, Blood vessel development 等富集

# Müller glia 上调基因 GO
muller_up_go <- do_go(muller_up_genes)
```

### 论文关键 GO 发现

- **RPE 上调**: 血管生成, VEGFA-VEGFR2 信号, 凋亡通路
- **RPE 下调**: 脂质储存, 离子稳态, 眼发育
- **Müller glia PR-like**: 光转导, 光敏感性

---

## Step 8: SASP 衰老评分

### 为什么做这一步？

AMD 是一种年龄相关疾病。SASP (Senescence-Associated Secretory Phenotype) 评估各细胞类型的**衰老程度**，判断衰老是否参与 AMD 发病。

```r
# 论文引用 Wang et al. [6] 的 SASP 基因列表
sasp_genes <- c("IL6", "IL8", "CXCL1", "CXCL2", "CXCL12", 
                "CCL2", "MMP1", "MMP3", "TIMP1", "TGFB1", "VEGFA", ...)

# 计算每个细胞类型的 SASP 评分
# 方法: 每个细胞中 SASP 基因的平均 log2(TPM) 表达值
DefaultAssay(integrated) <- "RNA"

# 获取 log2(TPM) 数据
expr <- GetAssayData(integrated, slot = "data")  # log-normalized data
sasp_expr <- expr[rownames(expr) %in% sasp_genes, ]

# 每个细胞的 SASP 评分 = 所有 SASP 基因的平均表达
integrated$sasp_score <- colMeans(sasp_expr, na.rm = TRUE)

# 按细胞类型和分组比较
VlnPlot(integrated, features = "sasp_score", group.by = "celltype", split.by = "group")

# 统计检验
for (ct in unique(integrated$celltype)) {
    sub <- subset(integrated, celltype == ct)
    amd_score <- sub$sasp_score[sub$group == "AMD"]
    ctrl_score <- sub$sasp_score[sub$group == "Control"]
    print(wilcox.test(amd_score, ctrl_score)$p.value)
}
```

---

## Step 9: 细胞间通讯分析 (CellPhoneDB)

### 为什么做这一步？

视网膜是一个多种细胞类型紧密互作的组织。CellPhoneDB 分析**配体-受体对**的表达，揭示 AMD 状态下细胞间通讯网络的重编程。

论文中使用的参数：
- 配体-受体对至少在一个细胞类型中 **10% 以上细胞**中表达
- 配体和受体均需可检测到表达
- 在细胞类型特异性层面比较 AMD vs Control 组的平均表达

```bash
# 安装 CellPhoneDB (Python)
pip install cellphonedb

# 准备输入文件
# 1. counts.txt: 基因表达矩阵 (normalized count)
# 2. meta.txt:   细胞注释文件 (barcode → celltype)
# meta.txt 格式:
#   Cell    cell_type
#   AAACCCA...    Rod
#   AAACGCT...    Cone

# 运行 CellPhoneDB (统计方法)
cellphonedb method statistical_analysis \
    meta.txt \
    counts.txt \
    --counts-data=gene_name \
    --threads=16 \
    --output-path=cellphonedb_out

# 可视化 (可选)
cellphonedb plot heatmap_plot meta.txt
cellphonedb plot dot_plot
```

### 输出解读

CellPhoneDB 输出每个配体-受体对在每对细胞类型之间的：
- **mean**: 配体和受体平均表达值的乘积
- **p-value**: 该互作是否显著高于随机背景

---

## Step 10: 可视化与结果汇报

### 关键可视化

```r
# 1. UMAP 细胞类型分布
DimPlot(integrated, group.by = "celltype", 
        cols = c("Rod"="#E41A1C", "Cone"="#377EB8", ...))

# 2. 细胞比例变化 (AMD vs Control)
cell_props <- integrated@meta.data %>%
    group_by(sample, group, celltype) %>%
    summarise(n = n()) %>%
    mutate(prop = n / sum(n))

ggplot(cell_props, aes(x = celltype, y = prop, fill = group)) +
    geom_boxplot() + theme_bw()

# 3. DEG 火山图
ggplot(sig_degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = celltype)) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed")

# 4. 关键基因表达小提琴图
VlnPlot(integrated, features = c("C3AR1", "KDR", "TF", "RBM3"), 
        group.by = "celltype", split.by = "group", pt.size = 0)

# 5. GO 富集气泡图
dotplot(rpe_up_go, showCategory = 10)

# 6. CellPhoneDB 互作网络热图
```

---

## 软件版本与参考文献

| 工具 | 版本 | 用途 |
|------|------|------|
| Cell Ranger | 6.x / 7.x | 比对与定量 |
| STAR | 2.7.x | RNA-seq 比对 (Cell Ranger 内置) |
| Seurat | 4.x / 5.x | 单细胞数据分析 |
| clusterProfiler | 4.x | GO 富集分析 |
| CellPhoneDB | 4.x / 5.x | 细胞间通讯 |
| R | 4.2+ | 统计计算环境 |
| Python | 3.9+ | CellPhoneDB 依赖 |

---

## 关键注意事项

1. **参考基因组**：食蟹猴 (Macaca fascicularis) 并非人类，需要下载对应的参考基因组和基因注释
2. **样本合并 vs 整合**：5 个样本来自 5 个独立 10X 文库，必须进行批次校正整合
3. **人工细胞注释**：自动注释工具 (如 SingleR) 可作为参考，但最终需要人工根据已知标记基因验证
4. **DEG 统计方法**：单细胞数据推荐 Wilcoxon 秩和检验 (非参数方法)，不假设正态分布
5. **多重检验校正**：所有 p 值需进行 Bonferroni 或 FDR 校正
6. **CellPhoneDB 需标准化表达值**：输入应为 normalized counts，非原始 UMI counts
