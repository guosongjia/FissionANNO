"""L3: ANNEVO ab initio prediction on residual regions + UniRef50 HGT filter.

Workflow:
  1. Hard-mask L1∪L2 gene regions (N)
  2. Extract residual intervals ≥ 1 kb
  3. ANNEVO prediction (Fungi model, CPU)
  4. ANNEVO decoding → GFF
  5. Remap coordinates to original genome
  6. L2 rescue: validate singleton_no_sog sidecar entries
  7. UniRef50 filter (HGT candidates only)
"""

rule l3_hard_mask:
    input:
        target = lambda w: SAMPLE_FASTA[w.sample],
        l1     = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.lifton.gff3"),
        l2     = os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.l2_kept.gff3"),
    output:
        masked = temp(os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.hardmasked.fasta")),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_hard_mask.log"),
    threads: 1
    shell:
        r"""
        python workflow/scripts/hard_mask_regions.py \
            --fasta {input.target} \
            --gffs {input.l1} {input.l2} \
            --output {output.masked} \
            > {log} 2>&1
        """

rule l3_extract_residuals:
    input:
        masked = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.hardmasked.fasta"),
    output:
        fasta = temp(os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.residuals.fasta")),
        maptsv = temp(os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.residuals.map.tsv")),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_extract_residuals.log"),
    params:
        min_bp = config["annevo"]["min_residual_bp"],
    threads: 1
    shell:
        r"""
        python workflow/scripts/extract_residuals.py \
            --fasta {input.masked} \
            --min-bp {params.min_bp} \
            --out-fasta {output.fasta} \
            --out-map {output.maptsv} \
            > {log} 2>&1
        """

rule l3_annevo_predict:
    input:
        fasta = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.residuals.fasta"),
    output:
        h5 = temp(os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.annevo_pred.h5")),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_annevo_predict.log"),
    params:
        annevo_dir = config["annevo"]["workflow_dir"],
        model      = config["annevo"]["model_path"],
        lineage    = config["annevo"]["lineage"],
        batch_size = config["annevo"]["batch_size"],
        env_bin    = config["annevo"]["conda_env_bin"],
    threads: config["resources"]["l3_annevo"]["threads"]
    shell:
        r"""
        export PATH="{params.env_bin}:$PATH"
        cd {params.annevo_dir}
        python prediction.py \
            -g {input.fasta} \
            -m {params.model} \
            -p {output.h5} \
            -l {params.lineage} \
            --batch_size {params.batch_size} \
            --num_workers {threads} \
            > {log} 2>&1
        """

rule l3_annevo_decode:
    input:
        fasta = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.residuals.fasta"),
        h5    = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.annevo_pred.h5"),
    output:
        gff = temp(os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.annevo_raw.gff")),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_annevo_decode.log"),
    params:
        annevo_dir = config["annevo"]["workflow_dir"],
        env_bin    = config["annevo"]["conda_env_bin"],
    threads: config["resources"]["l3_annevo"]["threads"]
    shell:
        r"""
        export PATH="{params.env_bin}:$PATH"
        cd {params.annevo_dir}
        python decoding.py \
            -g {input.fasta} \
            -p {input.h5} \
            -o {output.gff} \
            -t {threads} \
            > {log} 2>&1
        """

rule l3_remap:
    input:
        gff    = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.annevo_raw.gff"),
        maptsv = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.residuals.map.tsv"),
    output:
        gff = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.annevo_remapped.gff3"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_remap.log"),
    threads: 1
    shell:
        r"""
        python workflow/scripts/remap_coordinates.py \
            --gff {input.gff} \
            --map {input.maptsv} \
            --output {output.gff} \
            > {log} 2>&1
        """

rule l3_diamond:
    input:
        gff    = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.annevo_remapped.gff3"),
        genome = lambda w: SAMPLE_FASTA[w.sample],
    output:
        diamond = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.uniref50.tsv"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_diamond.log"),
    params:
        db     = config["uniref50"]["diamond_db"],
        evalue = config["uniref50"]["evalue"],
        min_prot = config["annevo"]["min_prot_len"],
    threads: config["resources"]["l3_uniref50"]["threads"]
    resources:
        mem_gb = config["resources"]["l3_uniref50"]["mem_gb"],
    shell:
        r"""
        python workflow/scripts/l3_diamond_search.py \
            --gff {input.gff} \
            --genome {input.genome} \
            --diamond-db {params.db} \
            --evalue {params.evalue} \
            --min-prot-len {params.min_prot} \
            --threads {threads} \
            --out-diamond {output.diamond} \
            > {log} 2>&1
        """

rule l3_rescue_l2:
    input:
        sidecar     = os.path.join(OUTDIR, "results", "{sample}", "sidecar", "{sample}.intra_genus_HGT_candidates.tsv"),
        miniprot    = os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.miniprot.gff3"),
        annevo_gff  = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.annevo_remapped.gff3"),
        diamond_tsv = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.uniref50.tsv"),
        l2_kept_in  = os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.l2_kept.gff3"),
        genome      = lambda w: SAMPLE_FASTA[w.sample],
    output:
        l2_kept_out = os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.l2_kept_rescued.gff3"),
        sidecar_out = os.path.join(OUTDIR, "results", "{sample}", "sidecar", "{sample}.sidecar_updated.tsv"),
        rescue_log  = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.rescue_log.tsv"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_rescue_l2.log"),
    params:
        min_aln_len = config["annevo"]["min_aln_len"],
    threads: 1
    shell:
        r"""
        python workflow/scripts/l2_rescue_from_l3.py \
            --sidecar {input.sidecar} \
            --miniprot-gff {input.miniprot} \
            --annevo-gff {input.annevo_gff} \
            --diamond-tsv {input.diamond_tsv} \
            --l2-kept-in {input.l2_kept_in} \
            --genome-fa {input.genome} \
            --l2-kept-out {output.l2_kept_out} \
            --sidecar-out {output.sidecar_out} \
            --rescue-log {output.rescue_log} \
            --min-aln-len {params.min_aln_len} \
            > {log} 2>&1
        """

rule l3_uniref50_filter:
    input:
        gff     = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.annevo_remapped.gff3"),
        genome  = lambda w: SAMPLE_FASTA[w.sample],
        l1      = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.lifton.gff3"),
        l2      = os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.l2_kept_rescued.gff3"),
        diamond = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.uniref50.tsv"),
    output:
        kept_gff = os.path.join(OUTDIR, "results", "{sample}", "L3", "{sample}.l3_kept.gff3"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l3_uniref50_filter.log"),
    params:
        keyword     = config["uniref50"]["schizo_taxon_keyword"],
        min_prot    = config["annevo"]["min_prot_len"],
        min_aln     = config["annevo"]["min_aln_len"],
        bam_dir     = config["annevo"]["bam_dir"],
        cov_floor   = config["annevo"]["coverage_floor"],
        schpo_cov   = config["annevo"]["schpo_min_target_cov"],
    threads: 1
    shell:
        r"""
        BAM="{params.bam_dir}/{wildcards.sample}.KAT_reads.bwa.sorted.bam"
        BAM_ARG=""
        if [ -f "$BAM" ]; then
            BAM_ARG="--bam $BAM --coverage-floor {params.cov_floor}"
        fi
        python workflow/scripts/l3_uniref50_filter.py \
            --braker-gff {input.gff} \
            --genome {input.genome} \
            --l1-gff {input.l1} \
            --l2-gff {input.l2} \
            --diamond-tsv {input.diamond} \
            --schizo-keyword "{params.keyword}" \
            --sample {wildcards.sample} \
            --threads {threads} \
            --min-prot-len {params.min_prot} \
            --min-aln-len {params.min_aln} \
            --schpo-min-target-cov {params.schpo_cov} \
            $BAM_ARG \
            --out-gff {output.kept_gff} \
            > {log} 2>&1
        """
