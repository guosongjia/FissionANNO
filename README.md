# FissionANNO

Population-scale annotation pipeline for *Schizosaccharomyces pombe* (extensible to other Schizosaccharomyces species).

See [CLAUDE.md](CLAUDE.md) for the design record from the grilling phase.

## Pipeline Overview

```mermaid
flowchart TD
    IN["<b>Input</b><br/>per-strain genome FASTA<br/>+ Lcon reference (FASTA + GFF3)"]

    subgraph L1["L1 — Reference transfer"]
        direction TB
        L1A["<b>lifton</b><br/>-sc 0.95  -copies<br/>-polish  -infer-genes"]
        L1B["<b>lifton_gff3_refine.py</b><br/>repair invalid ORFs:<br/>• missing stop → truncate to inframe stop (last 5%)<br/>&nbsp;&nbsp;&nbsp;or extend downstream<br/>• missing start → scan ≤ 300 bp upstream for ATG<br/>• un-refinable → pseudogene"]
        L1C["<b>build_unmapped_tsv.py</b><br/>4 reason classes:<br/>lifton_unmapped / refine_pseudogene /<br/>refine_truncated_at_contig_end /<br/>refine_frame_disrupted"]
        L1A --> L1B --> L1C
    end

    L1OUT[("L1/refine.gff3<br/>~5067–5103 genes/strain")]

    subgraph L2["L2 — Cross-species protein scan & SOG conflict resolution"]
        direction TB
        L2A["<b>miniprot</b> whole-genome scan<br/>query = 9 species protein FASTA<br/>(no masking)"]
        L2B["<b>① Pre-filter</b><br/>drop if any of:<br/>• aln &lt; 50 aa<br/>• identity &lt; 0.3<br/>• overlap with L1 gene ≥ 50%"]
        L2C["<b>② Locus collapse + SOG classification</b><br/>4 relation labels:<br/>• non_reference_gene<br/>• missing_lift<br/>• intra_genus_HGT_from_&lt;sp&gt;<br/>• diverged_paralog_or_misannot"]
        L2D["<b>③ HGT double gate</b> (HGT candidates only)<br/>candidate_id &gt; SOG Lcon×donor pairwise max<br/>(else demote → non_reference_gene)<br/>candidate_id ≥ 0.9 (else → sidecar)"]
        L2E["<b>④ ORF completeness filter</b><br/>require stop_codon<br/>CDS aa / query aa ≥ 0.95"]
        L2A --> L2B --> L2C --> L2D --> L2E
    end

    L2GFF[("l2_kept.gff3<br/>1–12 high-confidence<br/>new genes/strain")]
    L2TSV[("l2_candidates.tsv<br/>all 4 relation classes")]
    L2SC[("sidecar.tsv<br/>singleton_no_sog / low_conf_HGT /<br/>partial_no_stop / full_aln_no_stop /<br/>orf_too_short")]

    NEXT["L3 (BRAKER4 + UniRef50)<br/>or merge_layers → final.gff3"]

    IN --> L1
    L1 --> L1OUT
    L1OUT --> L2
    L2 --> L2GFF
    L2 --> L2TSV
    L2 --> L2SC
    L2GFF --> NEXT

    classDef io fill:#eef5ff,stroke:#3366aa,stroke-width:1px,color:#000
    classDef out fill:#e8f5e9,stroke:#2e7d32,stroke-width:1px,color:#000
    class IN,L1OUT io
    class L2GFF,L2TSV,L2SC out
```

**L1** transfers PomBase annotation gene-by-gene to the target strain, with
small-scale start/stop codon repair; unmapped genes are bucketed into 4 reason
classes. Main output: `refine.gff3` (~5k genes/strain).

**L2** uses a 9-species protein library to discover genes that L1 missed or
that were horizontally transferred within the genus. Only hits **not** overlapping
L1 are kept; the SOG table + pairwise identity assign one of 4 relation labels,
and an ORF-completeness filter pushes fragment noise to the sidecar. Main output:
1–12 high-confidence new genes per strain.

## Layout
```
FissionANNO/
  CLAUDE.md                    # decision record
  config/
    config.yaml                # all tunables
    manifest.tsv               # per-strain input
  workflow/
    Snakefile
    rules/
      common.smk
      l1_lifton.smk
      l2_miniprot.smk
      l3_braker4.smk
      merge.smk
    scripts/                   # python helpers
      lifton_gff3_refine.py    # in-tree copy
      build_sog_index.py
      build_unmapped_tsv.py
      l2_conflict_resolve.py
      softmask_regions.py
      run_braker4.py
      l3_uniref50_filter.py
      merge_layers.py
      capture_versions.py
    envs/
      lifton.yaml
      postprocess.yaml
  profiles/
    local/config.yaml          # 64-core single-machine profile
  resources/                   # cached intermediate (built once)
```

## Status
- 2026-05-22: scaffold + L1 refine fixes + SOG/unmapped scripts implemented
- 2026-05-23: conda env installed at `/data/c/jiaguosong/conda_envs/fissionanno` (899 MB)
- L2/L3/merge: rule wiring done; python scripts are placeholders (exit 2)

## Setup

Prerequisites: `micromamba` (or `mamba`) on PATH; `curl`, `tar`, `sed`.

```bash
git clone <repo> FissionANNO && cd FissionANNO
bash setup_env.sh                                                # default prefix
# or:
ENV_PREFIX=/path/to/your/envs/fissionanno bash setup_env.sh      # custom prefix
conda activate /data/c/jiaguosong/conda_envs/fissionanno
```

`setup_env.sh` is end-to-end reproducible: it neutralizes a stale
`~/.condarc`, builds the patched `cigar` wheel (lifton transitive dep
with broken upstream packaging), and installs lifton + pytest in the
right order. Tested from a clean state on 2026-05-23.

The repo is a single conda env. BRAKER4 is invoked as an *external*
Snakemake workflow at `/data/c/jiaguosong/BRAKER4/` via subprocess — it
is **not** installed into this env.

## Next
1. Run refine A1–A5 bug A/B test on 5 sample strains using existing `1.1_lifton_original_gff3` outputs.
2. Implement `build_sog_index.py`, `build_unmapped_tsv.py`, then `l2_conflict_resolve.py`.
3. Smoke-test L1 + L2 on the 5-strain manifest.
4. Wire `run_braker4.py` against the `/data/c/jiaguosong/BRAKER4` workflow.
5. Implement `l3_uniref50_filter.py` and `merge_layers.py`.

## Run
```bash
cd /data/c/jiaguosong/FissionANNO
snakemake --snakefile workflow/Snakefile --profile profiles/local -n   # dry run
```
