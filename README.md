# FissionANNO

A Snakemake pipeline for population-scale genome annotation in haploid fission yeasts (*Schizosaccharomyces*).

## Pipeline Overview

```mermaid
flowchart TD
    IN["<b>Input</b><br/>per-strain genome FASTA<br/>+ reference annotation (FASTA + GFF3)"]

    subgraph L1["L1 — Reference transfer"]
        direction TB
        L1A["<b>lifton</b><br/>-sc 0.95  -copies<br/>-polish  -infer-genes"]
        L1B["<b>lifton_gff3_refine.py</b><br/>repair invalid ORFs:<br/>• missing stop → truncate to inframe stop (last 5%)<br/>&nbsp;&nbsp;&nbsp;or extend downstream<br/>• missing start → scan ≤ 300 bp upstream for ATG<br/>• un-refinable → pseudogene"]
        L1C["<b>build_unmapped_tsv.py</b><br/>4 reason classes:<br/>lifton_unmapped / refine_pseudogene /<br/>refine_truncated_at_contig_end /<br/>refine_frame_disrupted"]
        L1A --> L1B --> L1C
    end

    L1OUT[("L1/refine.gff3")]

    subgraph L2["L2 — Cross-species protein scan & SOG conflict resolution"]
        direction TB
        L2A["<b>miniprot</b> whole-genome scan<br/>query = 9 species protein FASTA<br/>(no masking)"]
        L2B["<b>① Pre-filter</b><br/>drop if any of:<br/>• aln &lt; 50 aa<br/>• identity &lt; 0.3<br/>• overlap with L1 gene ≥ 50%"]
        L2C["<b>② Locus collapse + SOG classification</b><br/>4 relation labels:<br/>• non_reference_gene (absent from reference)<br/>• missing_lift (in reference but lifton failed)<br/>• intra_genus_HGT_from_&lt;sp&gt; (closer to donor than reference)<br/>• diverged_paralog_or_misannot (ambiguous)"]
        L2D["<b>③ HGT double gate</b> (HGT candidates only)<br/>candidate_id &gt; SOG ref×donor pairwise max<br/>(else demote → non_reference_gene)<br/>candidate_id ≥ 0.9 (else → sidecar)"]
        L2E["<b>④ ORF completeness filter</b><br/>require intact ORF<br/>CDS aa / query aa ≥ 0.95"]
        L2A --> L2B --> L2C --> L2D --> L2E
    end

    L2GFF[("l2_kept.gff3<br/>high-confidence new genes<br/>(non_ref / missing_lift / HGT)")]

    subgraph L3["L3 — Ab initio HGT discovery (ANNEVO)"]
        direction TB
        L3A["<b>hard_mask_regions.py</b><br/>N-mask all L1∪L2 gene regions"]
        L3B["<b>extract_residuals.py</b><br/>extract non-N intervals ≥ 1 kb"]
        L3C["<b>ANNEVO</b> (Fungi model, CPU)<br/>ab initio gene prediction with introns"]
        L3D["<b>UniRef50 annotation &amp; HGT filter</b><br/>DIAMOND blastp → exclude Schizo / TE / Metazoa<br/>contamination check: contig context + BAM coverage"]
        L3E["<b>L2 rescue</b><br/>validate singleton_no_sog → non_reference_gene"]
        L3A --> L3B --> L3C --> L3D --> L3E
    end

    L3GFF[("l3_kept.gff3<br/>putative extra-genus HGT<br/>+ L2 rescued genes")]

    subgraph MG["Merge"]
        direction TB
        MGA["<b>merge_layers.py</b><br/>L1 + L2 + L3 → sorted GFF3<br/>source= provenance tags"]
        MGB["<b>extract_proteins.py</b><br/>CDS → translated proteins<br/>gzipped FASTA"]
        MGA --> MGB
    end

    FINAL[("merged/final.gff3<br/>+ proteins.fa.gz")]

    IN --> L1
    L1 --> L1OUT
    L1OUT --> L2
    L2 --> L2GFF
    L2GFF --> L3
    L2E -.->|sidecar| L3E
    L3 --> L3GFF
    L1OUT --> MG
    L2GFF --> MG
    L3GFF --> MG
    MG --> FINAL

    classDef io fill:#eef5ff,stroke:#3366aa,stroke-width:1px,color:#000
    classDef out fill:#e8f5e9,stroke:#2e7d32,stroke-width:1px,color:#000
    class IN,L1OUT io
    class L2GFF,L3GFF,FINAL out
```

**L1** transfers the reference annotation gene-by-gene to the target strain using [LiftOn](https://github.com/Kuanhao-Chao/LiftOn), followed by a custom refinement step that repairs invalid ORFs (extending to a downstream stop codon, scanning upstream for a missing start codon, or reclassifying as pseudogene when repair is not possible). Genes that could not be transferred are catalogued with one of four reason classes. Main output: `refine.gff3`.

**L2** uses a multi-species protein library to discover genes that L1 missed or that were horizontally transferred within the genus. [miniprot](https://github.com/lh3/miniprot) scans the whole genome without masking; hits overlapping L1 genes are silently dropped. The remaining candidates are classified using a SOG (Schizosaccharomyces orthogroup) table: each locus is assigned one of four relation labels (*non_reference_gene*, *missing_lift*, *intra_genus_HGT_from_\<sp\>*, *diverged_paralog_or_misannot*). HGT candidates pass a double identity gate (absolute floor + SOG-derived pairwise cutoff), and an ORF-completeness filter pushes fragment noise to the sidecar (a separate TSV for manual review). Main output: 1–12 high-confidence new genes per strain.

**L3** targets extra-genus HGT in regions not covered by L1/L2. Gene-annotated regions are hard-masked (replaced with N), and residual intervals ≥ 1 kb are extracted and fed to [ANNEVO](https://github.com/xjtu-omics/ANNEVO), a pretrained deep-learning gene predictor with a Fungi-specific model that supports intron prediction without requiring retraining. Predictions are remapped to original genome coordinates and searched against UniRef50 with DIAMOND. A multi-layer filter removes hits to same-genus proteins, transposable elements, and likely assembly contaminants (verified by contig context and read-depth). As a side effect, L3 also rescues L2 sidecar entries (*singleton_no_sog*) that are independently validated by ANNEVO, reclassifying them as *non_reference_gene*. Main output: 0–3 putative extra-genus HGT candidates per strain.

**Merge** combines L1, L2, and L3 outputs into a single per-strain GFF3 sorted by genomic coordinate, with `source=` provenance tags on every gene and mRNA feature. A companion script translates all CDS features into a gzipped protein FASTA with provenance metadata in the headers, ready for downstream PAV/CNV/identity analyses.

## Setup

Prerequisites: `conda` on PATH.

```bash
git clone <repo> FissionANNO && cd FissionANNO
bash setup_env.sh
export PATH="/data/c/jiaguosong/conda_envs/fissionanno/bin:$PATH"
```

ANNEVO requires a separate environment (CPU-only setup):
```bash
conda create -p /data/c/jiaguosong/conda_envs/annevo python=3.10 -y
export PATH="/data/c/jiaguosong/conda_envs/annevo/bin:$PATH"
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install bcbio-gff h5py torchmetrics pandas numpy tqdm biopython
git clone https://github.com/xjtu-omics/ANNEVO.git /data/c/jiaguosong/ANNEVO
```

## Run
```bash
cd /data/c/jiaguosong/FissionANNO
export PATH="/data/c/jiaguosong/conda_envs/fissionanno/bin:$PATH"
snakemake --profile profiles/local -n          # dry run
snakemake --profile profiles/local --cores 64  # full run
```
