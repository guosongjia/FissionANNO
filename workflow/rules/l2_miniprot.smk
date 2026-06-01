"""L2: miniprot whole-genome scan with 9-species protein library, then SOG-based conflict resolution."""

rule l2_miniprot:
    input:
        target = lambda w: SAMPLE_FASTA[w.sample],
        proteins = os.path.join(OUTDIR, "shared", "combined_9sp.protein.fa"),
    output:
        gff = temp(os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.miniprot.gff3")),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l2_miniprot.log"),
    threads: lambda w: config["resources"]["l2_miniprot"]["threads"]
    resources:
        mem_gb = lambda w: config["resources"]["l2_miniprot"]["mem_gb"],
    params:
        extra = config["miniprot"]["extra_opts"],
    shell:
        r"""
        mkdir -p $(dirname {output.gff})
        miniprot -t {threads} --gff {params.extra} \
            {input.target} {input.proteins} \
            > {output.gff} 2> {log}
        """

rule l2_conflict_resolve:
    input:
        l1_gff   = os.path.join(OUTDIR, "results", "{sample}", "L1", "{sample}.refine.gff3"),
        l2_gff   = os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.miniprot.gff3"),
        sog_idx  = os.path.join(OUTDIR, "shared", "sog_index.pkl"),
        proteins = os.path.join(OUTDIR, "shared", "combined_9sp.protein.fa"),
        genome   = lambda w: SAMPLE_FASTA[w.sample],
    output:
        gff      = os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.l2_kept.gff3"),
        tsv      = os.path.join(OUTDIR, "results", "{sample}", "L2", "{sample}.l2_candidates.tsv"),
        sidecar  = os.path.join(OUTDIR, "results", "{sample}", "sidecar", "{sample}.intra_genus_HGT_candidates.tsv"),
    log:
        os.path.join(OUTDIR, "logs", "{sample}", "l2_resolve.log"),
    params:
        overlap_min   = config["conflict"]["overlap_reciprocal_min"],
        hgt_min_id    = config["conflict"]["hgt_min_identity"],
        diverged_max  = config["conflict"]["diverged_paralog_max_identity"],
        adjacent_max  = config["conflict"]["adjacent_max_distance_bp"],
        min_aln_aa    = config["conflict"]["min_aln_aa"],
        min_identity  = config["conflict"]["min_identity"],
        orf_min_cov   = config["conflict"]["orf_min_coverage"],
        ref_species   = config["sog"]["ref_species"],
    threads: lambda w: config["resources"]["l2_resolve"]["threads"]
    shell:
        r"""
        mkdir -p $(dirname {output.sidecar})
        python workflow/scripts/l2_conflict_resolve.py \
            --sample {wildcards.sample} \
            --l1-gff {input.l1_gff} \
            --l2-gff {input.l2_gff} \
            --sog-index {input.sog_idx} \
            --protein-fa {input.proteins} \
            --genome-fa {input.genome} \
            --overlap-min {params.overlap_min} \
            --hgt-min-identity {params.hgt_min_id} \
            --diverged-max-id {params.diverged_max} \
            --adjacent-max-bp {params.adjacent_max} \
            --min-aln-aa {params.min_aln_aa} \
            --min-identity {params.min_identity} \
            --orf-min-coverage {params.orf_min_cov} \
            --ref-species "{params.ref_species}" \
            --out-gff {output.gff} \
            --out-tsv {output.tsv} \
            --out-sidecar {output.sidecar} \
            > {log} 2>&1
        """
