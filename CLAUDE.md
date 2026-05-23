# FissionANNO — 项目决策记录

本文件归档 grilling 阶段（2026-05-20 ~ 2026-05-22）确定的所有架构与实现决策。后续讨论以本文件为基线，新决策追加到末尾的"变更日志"。

---

## 1. 项目目标与范围

- **目标**：为 ~830 株 *S. pombe* 自然分离株生成完整的每株基因组注释，包含 non-reference 基因（PomBase 漏注释、属内 HGT、候选属外 HGT）。
- **下游**：群体规模 PAV / CNV / identity_score 类分析，要求注释含 provenance 标签。
- **不做**：ncRNA/tRNA/rRNA（参考已剔除）、TE 单独建模、L4 ncRNA 层。

## 2. 三层架构（最终方案）

| 层 | 工具 | 输入 | 角色 |
|---|---|---|---|
| L1 | lifton | Lcon (= PomBase) | 参考骨架转移 |
| L2 | miniprot | 9-species 蛋白库（未去冗余）| 发现属内 HGT 与 PomBase 漏注释同源基因 |
| L3 | BRAKER4 (外部 subprocess) | soft-masked (L1∪L2) 残余区域 + UniRef50 后过滤 | 候选属外 HGT 与极少数残余基因 |

- L1 不做 strict mask，L2 全基因组扫描，按 SOG 表做冲突解决。
- L3 调用 `/data/c/jiaguosong/BRAKER4` 已有的 Snakemake 工作流，不重装。
- 候选基因没有 UniRef50 hit → 丢；top hit 非 Schizosaccharomyces → 标 `HGT_call=putative_<top_taxon>`。

## 3. 关键参数

### 3.1 lifton
- 去掉 `-exclude_partial`（让 partial 基因进 refine 自行处理）
- 保留 `-sc 0.95`（仅约束 `-copies` 副本，不影响主映射）
- 保留 `-copies`、`-polish`、`-infer-genes`
- 主映射阈值 `-s 0.5 -a 0.5` 维持 lifton 默认

### 3.2 refine 脚本
- 基本保留现有 `lifton_gff3_refine.py`，只做小尺度起止密码子修正
- 由数据决定的待修 bug（A1–A5）：用 5 株现有 1.1 输出做修复前后对照测试，差异判优后再合并
  - A1: `extend_to_upstream_start` 短窗口分支返回值数量不匹配
  - A2: `setup_logging` 不挂 handler，warning/info 全丢
  - A3: `match_exon_cds` 修改副作用 + 返回 False 引发双重重建
  - A4: `extend_to_downstream_stop` 边界 codon 仅检 stop，破坏 reading frame
  - A5: `extend_to_upstream_start` 未保证 `needed % 3 == 0`
- 参数化（B1/B3）：`--last-pct-frac`（默认 0.05）+ 绝对下限；`--upstream-window-bp` 默认 `min(L/2, 300)`
- B2 跨 intron 延伸限制保留但写入 stat 表
- B4 pseudogene 复活逻辑保留现状（不修）
- 输出（C1–C3）：stat 加表头；GFF 加 `##gff-version 3`；流程结尾对所有 CDS 强制 phase recompute

### 3.3 L2 miniprot
- 蛋白库 = 9 物种 `*.protein.fasta` 直接 cat（option A，保留 in-paralog 信号），各 ID 加物种前缀
- 全基因组扫描，无 strict mask
- singleton / SOG 缺失蛋白同样跑，但 hit 进 sidecar `intra_genus_HGT_candidates.tsv`，不进主 GFF

### 3.4 L3 BRAKER4
- 外部 subprocess 调用 `/data/c/jiaguosong/BRAKER4/run_snakemake.sh` 风格
- 蛋白证据 = `/data/c/jiaguosong/BRAKER4/orthodb/Fungi.9FissionYeast_251010.fa`
- 无 RNA-seq（830 株均为天然分离株，无种群级 RNA-seq）
- 后过滤：UniRef50 DIAMOND，**有 hit 即留**；top hit 非 Schizosaccharomyces → `HGT_call=putative`

## 4. 冲突解决规则（L1 / L2）

### Case 1：M、L 空间重叠（≥50% reciprocal）
- `OG_M == OG_L` → 丢 M
- 同 family / 嵌套 SOG → 标 `paralog_suspect`，保留双方
- 完全无关 → 标 `real_conflict`

### Case 2：M 不与任何 L 重叠
计算 P_M 对 P_lcon、P_other_i 的 identity：
- 最近邻 = P_lcon → `missing_lift`
- 最近邻 = 某 P_other_j 且 identity 优势 ≥ 10pp（或 bit-score ≥ 50） → `intra_genus_HGT_from_<species_j>`
- 全部 < 60% → `diverged_paralog_or_misannot`，hold for review
- Lcon 在 OG_M 无成员 → `non_reference_gene`（金标准发现）

### Case 3：邻近 < 2 kb 不重叠（坍缩）
- 都完整 → 当 Case 2 处理
- L 完整 / M 片段 → 丢 M
- M 完整 / L 片段 → 保留双方，L 标 `lift_partial_superseded_by_M`
- 都片段 → `ambiguous_fragment_pair`

完整性判据：refine 输出的 `valid_orf` 属性。

阈值 10pp / 50bit 为暂定，后续按 ROC 调。

## 5. 资源与环境

- 服务器：单机 64 核 / 251 GB / `/data/c` 905 GB（`/data/a` 已满，禁止写入）
- 调度：Snakemake `--cores 64`，per-rule threads 见配置
- conda 环境拆分：`envs/lifton.yaml`（lifton + miniprot + gffutils + biopython）、`envs/postprocess.yaml`（pyranges + pandas + diamond）、`envs/braker_caller.yaml`（仅 apptainer 客户端）
- BRAKER4 不进 conda；走外部 subprocess
- 仓库路径：`/data/c/jiaguosong/FissionANNO/`
- 默认 outdir：`/data/c/jiaguosong/ngs_assembly/Spombe/annotation`（config 可覆盖）

## 6. 资源文件路径

- SOG 表：`/data/a/jiaguosong/jiaguosong/reference_genome_panel/nine_representative_strains_related/260511_SOG_updates_and_related_analysis/1_SOG_table_update/SOG_table.260511.updated2.tsv`
- 9 物种蛋白 FASTA：`/data/a/jiaguosong/jiaguosong/FYRP/annotation/*.protein.fasta`
- UniRef50 DIAMOND：`/data/c/resource/database/uniref50_2022_12_14.dmnd`
- BRAKER4 工作流：`/data/c/jiaguosong/BRAKER4/`
- BRAKER4 蛋白库：`/data/c/jiaguosong/BRAKER4/orthodb/Fungi.9FissionYeast_251010.fa`

参考基因组路径**必须 config 化**——用户后续会替换 reference annotation。

## 7. 输入与输出

### 7.1 Manifest（3 列 TSV）
```
sample_name<TAB>genome_fasta<TAB>species
```
`species` 默认 `Schizosaccharomyces_pombe`，用于多物种扩展。当前 manifest 路径 `/data/c/jiaguosong/ngs_assembly/Spombe/lifton_transfer/pombe_strain` 是 2 列，迁移时补全 species 列。

### 7.2 输出布局（Snakemake 风格）
```
{outdir}/
  results/{sample}/
    L1/{sample}.lifton.gff3
    L1/{sample}.refine.gff3
    L1/{sample}.stat.tsv
    L2/{sample}.miniprot.gff3
    L2/{sample}.conflict_resolution.tsv
    L3/{sample}.softmasked.fasta
    L3/{sample}.braker4.gff3
    L3/{sample}.uniref50.tsv
    merged/{sample}.final.gff3
    sidecar/{sample}.intra_genus_HGT_candidates.tsv
    logs/{sample}.{rule}.log
  unmapped_summary.tsv          # all-strain aggregate
  versions.tsv                  # tool versions per run
```

### 7.3 unmapped_summary schema
```
sample_name<TAB>gene_id<TAB>reason<TAB>parent_seq<TAB>note
```
reason 取值（4 类，由 `build_unmapped_tsv.py` 实际生成）：
- `lifton_unmapped` — lifton 完全没找到映射（来自 `stats/unmapped_features.txt`）
- `refine_pseudogene` — refine 标 pseudogene
- `refine_truncated_at_contig_end` — refine 救不回 + 在 contig 边界（assembly 不完整，~99% 此类）
- `refine_norf_uncorrectable` — refine 救不回 + 非边界（真 ORF 损伤，~1% 此类）

note 列若有 lifton flag 信息（如 `frameshift;stop_codon_gain`）则填入。

### 7.4 Provenance 强制标签
- `source=lifton|miniprot_L2|braker4_L3`
- `evidence=UniRef50_<acc>` / `SOG_id=<id>` / `HGT_call=...`
- `relation=missing_lift|HGT_candidate|paralog_suspect|...`

## 8. 测试集（5 株）
DY47073、DY46687、DY44518-zxr、DY42495-zxr、DY39827。开发期所有 rule 在测试集跑通后再上全量 830 株。

## 9. 工程要求

- 每株独立 workdir，`ln -s` 替代 cp；跑完移最终产物到 `results/{sample}/`，删除 workdir
- 错误隔离：单株失败不阻断其他株
- 断点续跑：Snakemake 原生自动跳过已存在 output
- 每株 log 记录工具版本、完整命令行、输入 fasta md5、起止时间戳，聚合到 `versions.tsv`
- 参数全部走 `config/config.yaml`，关键参数命令行可覆盖
- 参考基因组路径**绝不硬编码**

## 10. 已知未决与延后项

- HGT identity-advantage 阈值 10pp / 50bit：依经验，后续 ROC 调
- L3 候选属外 HGT 真假判定：仅做 UniRef50 标注，最终判断留人工
- 用户后续会提供新的 reference genome annotation：架构需对参考切换零代码改动
- A1–A5 修复决策依赖 5 株对照测试结果
- **SOG 表覆盖度**：Lcon 5186 个基因中 5027 个（97%）在 SOG 表内；**159 个 Lcon 基因不在任何 SOG**。L2 冲突解决对这部分基因采用 fallback (b')：**无条件保留 M，标 `relation=no_sog_info_lift_unverifiable`**，不做阈值判定，由下游人工裁决。
- **SOG 表数据 bug**：source TSV line 1634 (SOG_5834) 列数 10 vs 期望 13（少 3 个 tab）。已于 2026-05-22 修复，备份在 `*.bak.20260522`。`build_sog_index.py` 默认 strict 模式会 fail loud，`--no-strict` 可跳过坏行。
- **MT 基因预过滤**：reference GFF 中的 `SPMIT.*`（8 个 mitochondrial genes）会被 lifton 错误转到核 contig（assembly 已去除 MT），产生假救回（如 SPMIT.11 → 10 aa 假 ORF）。**用户负责**在交付新 reference GFF 时预过滤掉 MT 基因；FissionANNO 不再单独处理。
- **去 `-exclude_partial` 的真实收益（5 株统计）**：每株多 23-34 个基因；其中 ~5% 真实救回（如 SPBCPT2R1.10 经下游延伸 132 bp 找到新 stop，与 ref 99.1% identity）、~26% 正确归类为 pseudogene、~70% 标 truncated_at_contig_end（assembly 截断）。整体修改方向正确。

## 变更日志
- 2026-05-22：initial commit，记录 grilling 第 1 轮所有决策。
- 2026-05-22：refine 脚本测试与修复完成；A1+A2 修（确认 5 株上零差异）、B1 (last_pct_min_bp=30) 救回 0–3 个短伪基因/株、B3 (upstream_window_bp=300) 5 株上未触发但作为防御保留、C1/C2/C3 加表头与 phase 兜底。`--no-final-phase-recompute` 提供 legacy 兼容开关。
- 2026-05-22：实施 `build_sog_index.py`，发现 SOG 表 line 1634 列数错误（已记入第 10 节）；Lcon 5027/5186 基因有 SOG 映射，159 基因需 fallback 规则。
- 2026-05-22：修复 SOG 表 line 1634 数据 bug（备份 `*.bak.20260522`），SOG 总数 5311 → 5312。strict 模式现可通过。
- 2026-05-23：实施 `build_unmapped_tsv.py`；refine stat 表加 `truncated_orf` 列；unmapped_summary 简化为 4 类 reason（`lifton_unmapped` / `refine_pseudogene` / `refine_truncated_at_contig_end` / `refine_norf_uncorrectable`）。`lifton_low_score`(damaging) 不作为独立 reason，由 score.txt 单独输出非健康部分。5 株分析显示 ~99% 的 norf_no_pseu 是 contig-end truncation（assembly 截断），仅 ~1% 是真受损基因。
- 2026-05-23：env 安装完成（micromamba + 单环境，899 MB），加 `setup_env.sh` 端到端可重现脚本，git init。
- 2026-05-23：L1 端到端 5 株跑通；修复两个 rule bug（lifton -g 改 ref_db 避免 gff3_db 并发竞争；unmapped 路径加 L1/）。每株 5067-5103 个基因，符合预期。识别两个待办：MT 基因 lifton 假救回（用户预过滤）、SPBCPT2R1.10 类自然变异下游延伸救回（refine 行为正确，保留）。
