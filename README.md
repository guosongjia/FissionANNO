# FissionANNO

Population-scale annotation pipeline for *Schizosaccharomyces pombe* (extensible to other Schizosaccharomyces species).

See [CLAUDE.md](CLAUDE.md) for the design record from the grilling phase.

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

```bash
# 1. Patch cigar source (lifton transitive dep; upstream packaging is broken)
cd /tmp && curl -sL https://files.pythonhosted.org/packages/source/c/cigar/cigar-0.1.3.tar.gz | tar xz
cd cigar-0.1.3 && sed -i '/^import ez_setup/d; /^ez_setup.use_setuptools/d' setup.py
/data/c/jiaguosong/conda_envs/fissionanno/bin/pip install --no-build-isolation --no-deps .

# 2. Create env (skips cigar; rest of pip section installs)
CONDARC=/tmp/fa_condarc_dir/condarc \
  micromamba env create -y --override-channels \
    -c conda-forge -c bioconda \
    -p /data/c/jiaguosong/conda_envs/fissionanno \
    -f workflow/envs/fissionanno.yaml

# 3. Activate
conda activate /data/c/jiaguosong/conda_envs/fissionanno
```

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
