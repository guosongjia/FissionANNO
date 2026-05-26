"""Merge L1 + L2_kept_rescued + L3_kept into a single per-strain GFF with provenance tags."""

rule merge_per_strain:
    input:
        l1 = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.refine.gff3"),
        l2 = os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.l2_kept_rescued.gff3"),
        l3 = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.l3_kept.gff3"),
    output:
        gff = os.path.join(OUTDIR, "results", "{sample}", "merged", "{sample}.final.gff3"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "merge.log"),
    threads: lambda w: config["resources"]["merge"]["threads"]
    shell:
        r"""
        python workflow/scripts/merge_layers.py \
            --sample {wildcards.sample} \
            --l1 {input.l1} \
            --l2 {input.l2} \
            --l3 {input.l3} \
            --output {output.gff} \
            > {log} 2>&1
        """
