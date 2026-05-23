"""L1: lifton transfer + refine."""

rule l1_lifton:
    input:
        target  = lambda w: SAMPLE_FASTA[w.sample],
        ref_fa  = config["reference"]["fasta"],
        ref_gff = config["reference"]["gff3"],
        ref_db  = config["reference"]["gff3_db"],
        prot    = config["reference"]["protein_fa"],
        cds     = config["reference"]["cds_fa"],
    output:
        gff = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.lifton.gff3"),
        outdir_tag = os.path.join(OUTDIR, "results", "{sample}", "L1", ".lifton_done"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l1_lifton.log"),
    params:
        sc = config["lifton"]["sc"],
        s  = config["lifton"]["s"],
        a  = config["lifton"]["a"],
        flags = lambda w: " ".join([
            "-copies"        if config["lifton"]["copies"] else "",
            "-polish"        if config["lifton"]["polish"] else "",
            "-infer-genes"   if config["lifton"]["infer_genes"] else "",
            "-exclude_partial" if config["lifton"]["exclude_partial"] else "",
        ]).strip(),
        workdir = lambda w: os.path.join(OUTDIR, "results", w.sample, "L1", "_work"),
    threads: lambda w: config["resources"]["l1_lifton"]["threads"]
    resources:
        mem_gb = lambda w: config["resources"]["l1_lifton"]["mem_gb"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.workdir}
        ln -sf {input.target} {params.workdir}/{wildcards.sample}.fasta
        ln -sf {input.ref_fa} {params.workdir}/$(basename {input.ref_fa})
        cd {params.workdir}
        lifton \
            -g {input.ref_gff} \
            -P {input.prot} \
            -T {input.cds} \
            -o {wildcards.sample}.lifton.gff3 \
            -s {params.s} -a {params.a} -sc {params.sc} \
            {params.flags} \
            -t {threads} \
            {wildcards.sample}.fasta $(basename {input.ref_fa}) \
            > {log} 2>&1
        mv {wildcards.sample}.lifton.gff3 {output.gff}
        # Preserve lifton_output for unmapped extraction
        if [ -d lifton_output ]; then
            mv lifton_output ../{wildcards.sample}.lifton_output
        fi
        touch {output.outdir_tag}
        """

rule l1_refine:
    input:
        gff = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.lifton.gff3"),
        target = lambda w: SAMPLE_FASTA[w.sample],
    output:
        gff  = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.refine.gff3"),
        stat = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.stat.tsv"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l1_refine.log"),
    params:
        last_pct_frac      = config["refine"]["last_pct_frac"],
        last_pct_min_bp    = config["refine"]["last_pct_min_bp"],
        upstream_window_bp = config["refine"]["upstream_window_bp"],
        log_level          = config["refine"]["log_level"],
    threads: 1
    shell:
        r"""
        python workflow/scripts/lifton_gff3_refine.py \
            -g {input.gff} \
            -f {input.target} \
            -o {output.gff} \
            --log-file {output.stat} \
            --log-level {params.log_level} \
            > {log} 2>&1
        """

rule l1_unmapped_per_sample:
    input:
        lifton_outdir = os.path.join(OUTDIR, "results", "{sample}", "L1", ".lifton_done"),
        stat = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.stat.tsv"),
        gff  = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.refine.gff3"),
    output:
        tsv = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.unmapped.tsv"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l1_unmapped.log"),
    shell:
        r"""
        python workflow/scripts/build_unmapped_tsv.py \
            --sample {wildcards.sample} \
            --lifton-outdir {OUTDIR}/results/{wildcards.sample}/{wildcards.sample}.lifton_output \
            --refine-stat {input.stat} \
            --refined-gff {input.gff} \
            --output {output.tsv} \
            > {log} 2>&1
        """

rule unmapped_summary:
    input:
        per_sample = expand(os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.unmapped.tsv"),
                            sample=SAMPLES),
    output:
        tsv = os.path.join(OUTDIR, "unmapped_summary.tsv"),
    shell:
        r"""
        echo -e "sample_name\tgene_id\treason\tparent_seq\tnote" > {output.tsv}
        for f in {input.per_sample}; do tail -n +2 "$f" >> {output.tsv}; done
        """
