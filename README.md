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
- 2026-05-22: scaffold only. All `workflow/scripts/*.py` except `lifton_gff3_refine.py` are placeholders that exit 2.
- L1 rules wired against in-tree refine script. L2/L3/merge rules wired but call stubs.

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
