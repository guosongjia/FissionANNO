"""Shared helpers: protein library construction, version logging, manifest checks."""

# Combined 9-species protein library (built once, reused for all samples)
rule l2_build_combined_proteins:
    output:
        fa = os.path.join(OUTDIR, "shared", "combined_9sp.protein.fa"),
    params:
        lib_dir = config["miniprot"]["protein_lib_dir"],
        files   = config["miniprot"]["species_files"],
    log:
        os.path.join(OUTDIR, "logs", "shared", "build_combined_proteins.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.fa})
        : > {output.fa}
        for f in {params.files}; do
            sp=${{f%.protein.fasta}}
            awk -v sp="$sp" '/^>/{{sub(/^>/, ">" sp "|"); print; next}} {{print}}' \
                {params.lib_dir}/$f >> {output.fa}
        done 2> {log}
        """

# SOG table -> protein_id => SOG_id mapping (pickle)
rule l2_build_sog_index:
    input:
        sog = config["sog"]["path"],
        proteins = os.path.join(OUTDIR, "shared", "combined_9sp.protein.fa"),
    output:
        idx = os.path.join(OUTDIR, "shared", "sog_index.pkl"),
    params:
        ref_species = config["sog"]["ref_species"],
    log:
        os.path.join(OUTDIR, "logs", "shared", "build_sog_index.log"),
    shell:
        r"""
        python workflow/scripts/build_sog_index.py \
            --sog {input.sog} \
            --protein-fa {input.proteins} \
            --ref-species "{params.ref_species}" \
            --output {output.idx} \
            2> {log}
        """

# Version capture (run once, after `all`)
rule capture_versions:
    output:
        tsv = os.path.join(OUTDIR, "versions.tsv"),
    shell:
        r"""
        python workflow/scripts/capture_versions.py > {output.tsv}
        """
