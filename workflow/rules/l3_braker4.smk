"""L3: soft-mask L1∪L2 regions, run external BRAKER4 workflow, filter by UniRef50.

BRAKER4 is invoked as an external Snakemake subprocess against /data/c/jiaguosong/BRAKER4
(see CLAUDE.md §3.4). Do not attempt to install BRAKER4 inside this conda env.
"""

rule l3_softmask:
    input:
        target = lambda w: SAMPLE_FASTA[w.sample],
        l1     = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.refine.gff3"),
        l2     = os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.l2_kept.gff3"),
    output:
        masked = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.softmasked.fasta"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_softmask.log"),
    threads: lambda w: config["resources"]["l3_mask"]["threads"]
    shell:
        r"""
        python workflow/scripts/softmask_regions.py \
            --fasta {input.target} \
            --gffs {input.l1} {input.l2} \
            --output {output.masked} \
            > {log} 2>&1
        """

rule l3_braker4:
    input:
        masked   = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.softmasked.fasta"),
        protein  = config["braker4"]["protein_fa"],
    output:
        gff = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.braker4.gff3"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_braker4.log"),
    params:
        wf_dir   = config["braker4"]["workflow_dir"],
        busco    = config["braker4"]["busco_lineage"],
        enabled  = config["braker4"]["enabled"],
    threads: lambda w: config["resources"]["l3_braker4"]["threads"]
    resources:
        mem_gb = lambda w: config["resources"]["l3_braker4"]["mem_gb"],
    shell:
        r"""
        if [ "{params.enabled}" != "True" ]; then
            mkdir -p $(dirname {output.gff})
            : > {output.gff}
            echo "BRAKER4 disabled via config" > {log}
            exit 0
        fi
        python workflow/scripts/run_braker4.py \
            --sample {wildcards.sample} \
            --masked-fasta {input.masked} \
            --protein {input.protein} \
            --busco-lineage {params.busco} \
            --braker4-workflow {params.wf_dir} \
            --threads {threads} \
            --output-gff {output.gff} \
            > {log} 2>&1
        """

rule l3_uniref50_filter:
    input:
        gff = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.braker4.gff3"),
        masked = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.softmasked.fasta"),
    output:
        kept_gff = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.l3_kept.gff3"),
        diamond  = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.uniref50.tsv"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_uniref50.log"),
    params:
        db      = config["uniref50"]["diamond_db"],
        evalue  = config["uniref50"]["evalue"],
        keyword = config["uniref50"]["schizo_taxon_keyword"],
    threads: lambda w: config["resources"]["l3_uniref50"]["threads"]
    resources:
        mem_gb = lambda w: config["resources"]["l3_uniref50"]["mem_gb"],
    shell:
        r"""
        python workflow/scripts/l3_uniref50_filter.py \
            --braker-gff {input.gff} \
            --genome {input.masked} \
            --diamond-db {params.db} \
            --evalue {params.evalue} \
            --schizo-keyword "{params.keyword}" \
            --threads {threads} \
            --out-gff {output.kept_gff} \
            --out-diamond {output.diamond} \
            > {log} 2>&1
        """
