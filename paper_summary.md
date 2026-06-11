# Single-cell transcriptome atlas of spontaneous dry age-related macular degeneration in macaques

**期刊**: Fundamental Research 5 (2025) 1034–1046  
**DOI**: [10.1016/j.fmre.2023.02.028](https://doi.org/10.1016/j.fmre.2023.02.028)  
**通讯作者**: Wenru Su, 中山大学中山眼科中心, 眼科学国家重点实验室

---

## 一、研究背景

- **年龄相关性黄斑变性 (AMD)** 是老年人致盲的主要原因之一
- AMD 分为**干性 (dry/non-neovascular)** 和**湿性 (wet/neovascular)** 两种形式
- 干性 AMD 占 AMD 病例的 ~85-90%，目前缺乏有效的治疗方法
- **非人灵长类 (NHP)** 动物模型与人类视网膜高度相似，是研究 AMD 的理想模型
- 单细胞 RNA 测序 (scRNA-seq) 技术为在单细胞分辨率下解析疾病机制提供了有力工具

---

## 二、研究目的

利用 scRNA-seq 技术构建**自发性干性 AMD 食蟹猴**视网膜的单细胞转录组图谱，揭示干性 AMD 的细胞类型特异性转录组改变及其发病机制。

---

## 三、实验方法

### 3.1 动物模型筛选
- 对 **>600 只**食蟹猴进行眼科筛查（眼压、裂隙灯、眼底照相、OCT）
- 纳入标准：黄斑区出现玻璃膜疣 (drusen)，OCT 确认位于 RPE 层
- 最终入选 **6 只**自发性干性 AMD 猴 + **4 只**健康老龄猴（均为雌性, 19-22 岁）

### 3.2 单细胞测序与数据分析
- 10X Genomics 单细胞平台 (Single Cell 3' Reagent Kits v3)
- 每样本目标捕获 10,000 个细胞
- 质控后保留约 **30,000** 个高质量视网膜细胞
- 鉴定出 **9 种细胞类型**：

| 细胞类型 | 英文名 | 所属层次 |
|---------|--------|---------|
| 视杆细胞 | Rod cells | 神经视网膜 |
| 视锥细胞 | Cone cells | 神经视网膜 |
| 双极细胞 | Bipolar cells (BCs) | 神经视网膜 |
| 无长突细胞 | Amacrine cells (ACs) | 神经视网膜 |
| 水平细胞 | Horizontal cells (HCs) | 神经视网膜 |
| 视网膜神经节细胞 | Retinal ganglion cells (RGCs) | 神经视网膜 |
| Müller 胶质细胞 | Müller glial cells | 神经视网膜 |
| 小胶质细胞 | Microglia | 神经视网膜 |
| 视网膜色素上皮细胞 | RPE cells | 脉络膜层 |

### 3.3 其他分析方法
- **CellPhoneDB**: 细胞间通讯分析
- **SASP 评分**: 衰老相关分泌表型评分
- **GO 富集分析**: 功能通路分析
- **小鼠模型验证**: NaIO₃ 诱导干性 AMD 小鼠模型 + Fer-1 (铁死亡抑制剂) 治疗
- **ERG**: 视网膜电图功能检测
- **H&E 染色**: 组织学分析

---

## 四、主要结果

### 4.1 单细胞图谱的建立
- 成功构建正常老龄与 SD-AMD 食蟹猴的视网膜单细胞转录组图谱
- SD-AMD 组中**视锥、视杆和 RPE 细胞比例下降**，提示这些细胞类型参与 AMD 发病

### 4.2 全局转录组改变
- BCs、视杆细胞和 Müller 胶质细胞中**上调 DEGs 较多**
- HCs、RGCs 和 BCs 中**下调 DEGs 较多**
- 提示这些细胞类型受 AMD 影响最为显著

### 4.3 Müller 胶质细胞的两个亚型
- 发现一个**类光感受器 (PR-like)** Müller 胶质亚型
  - 表达光感受器标志基因: **RHO, GNAT1, ROM1**
  - 可能具有向光感受器分化的潜能
- SD-AMD 组中该亚型比例**显著降低**
- AMD Müller 胶质中 **KDR (VEGFR2)** 表达上调，提示新生血管相关改变

### 4.4 小胶质细胞的激活
- SD-AMD 小胶质细胞表现出**激活特征**
- **C3AR1** 表达上调，与补体系统激活相关
- 观察到 **M2 型极化**标志物上调（促血管生成表型）
- 提示小胶质细胞参与干性 AMD 的免疫病理过程

### 4.5 光感受器的差异性改变

**视杆细胞**:
- 光敏感性和功能相关基因上调
- 表现出对氧化应激的代偿性反应

**视锥细胞**:
- 细胞周期相关基因 **GINS2, RTN** 下调
- **ABCC4** (cAMP 转运相关) 下调
- 细胞类型特异性转录组改变显著

### 4.6 RPE 细胞的新生血管特征

**上调基因**:
- **RAC1** — 内皮细胞侵袭调控
- **PTN** — 多效生长因子
- **CLIC4** — 氯离子通道
- **EGFL7** — EGF 样结构域蛋白
- **TGFBR3** — TGF-β 受体
- **STRA6** — 维生素 A 转运关键调控因子

**上调通路**:
- 血管生成 (Angiogenesis)
- 血管发育 (Blood vessel development)
- VEGFA-VEGFR2 信号通路
- 凋亡相关通路

**下调通路**:
- 脂质储存与定位
- 离子稳态与阳离子稳态
- 眼发育相关通路

> **重要发现**: 即使在早期干性 AMD 阶段，RPE 已表现出新生血管相关的转录组改变，这可能决定疾病向湿性 AMD 的转化。

### 4.7 铁死亡 (Ferroptosis) 在 AMD 中的作用

**四类神经细胞 (RGCs, HCs, ACs, BCs) 的共同改变**:
- 共同上调基因: **RBM3, CIRBP** (冷休克蛋白/应激反应)
- 共同下调基因: **NCDN, TF (转铁蛋白), APLP1**
- **转铁蛋白 TF 下调** → 铁代谢紊乱 → 铁积累 → **铁死亡**

**小鼠模型验证**:
- NaIO₃ 诱导干性 AMD 小鼠
- **Fer-1 (铁死亡抑制剂)** 治疗后显示出保护效果

---

## 五、讨论与结论

### 5.1 核心发现总结
1. **首次**在非人灵长类模型中构建了干性 AMD 的单细胞转录组图谱
2. Müller 胶质细胞的 PR-like 亚型可能是**视网膜再生的潜在细胞来源**，该亚型在 AMD 中减少
3. 小胶质细胞通过 **C3AR1** 介导的补体激活和 M2 极化参与 AMD 发病
4. RPE 在干性 AMD 早期即表现出**新生血管相关转录组改变**
5. 神经细胞中**铁死亡通路**的激活是干性 AMD 的共同特征

### 5.2 研究意义
- 为理解干性 AMD 的发病机制提供了**单细胞分辨率的路线图**
- 鉴定了多个**潜在治疗靶点**: C3AR1, KDR/VEGFR2, 铁死亡通路
- Müller 胶质细胞的再生潜力为**细胞治疗**提供了新思路

### 5.3 局限性
- 样本量相对较小 (6 AMD + 4 对照)
- 仅纳入了雌性食蟹猴
- 需要更大样本和功能性验证实验

---

## 六、关键基因与通路汇总

| 类别 | 基因/分子 | 变化方向 | 功能意义 |
|------|----------|---------|---------|
| Müller 胶质标记 | RHO, GNAT1, ROM1 | PR-like 亚型中表达 | 光感受器特征 |
| Müller 胶质 AMD | KDR (VEGFR2) | ↑ | 新生血管信号 |
| 小胶质细胞 | C3AR1 | ↑ | 补体系统激活 |
| RPE 新生血管 | RAC1, PTN, CLIC4, EGFL7 | ↑ | 血管生成 |
| RPE 代谢 | STRA6 | ↑ | 维生素 A 转运 |
| 视锥细胞 | GINS2, ABCC4 | ↓ | 功能障碍 |
| 铁死亡(共同) | RBM3, CIRBP | ↑ | 应激反应 |
| 铁死亡(共同) | NCDN, TF, APLP1 | ↓ | 铁代谢紊乱 |

---

## 参考文献 (精选)

1. Jonas JB, et al. Updates on the epidemiology of age-related macular degeneration. *Asia Pac J Ophthalmol*. 2017.
2. Mitchell P, et al. Age-related macular degeneration. *Lancet*. 2018.
3. Wang S, et al. Deciphering primate retinal aging at single-cell resolution. (引用 #7-8)
4. Pollak J, et al. ASCL1 reprograms mouse Müller glia into neurogenic retinal progenitors. *Development*. 2013.
5. Hahn P, et al. Maculas affected by AMD contain increased chelatable iron. *Arch Ophthalmol*. 2003.

---

> **论文全文**: 共 55 篇参考文献，发表于 Fundamental Research 2025 年第 5 卷第 3 期，页码 1034-1046。
