rule build_external_sample_set:
    input:
        fam=rules.standardize_genotype.output.fam,
    output:
        keep=ws_path(
            "external_samples/eligible_samples.keep"
        ),
        exclusions=ws_path(
            "external_samples/external_sample_exclusions.tsv"
        ),
        report=ws_path(
            "external_samples/external_samples.report.tsv"
        ),
    conda:
        "../envs/r_environment.yaml"
    params:
        keep_files=lambda wc: ";".join(
            map(
                str,
                cfg("external_samples/keep_files", []),
            )
        ),
        remove_files=lambda wc: ";".join(
            map(
                str,
                cfg("external_samples/remove_files", []),
            )
        ),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.keep}")"

        Rscript workflow/scripts/build_external_sample_set.R \
            --fam "{input.fam}" \
            --keep-output "{output.keep}" \
            --exclusions-output "{output.exclusions}" \
            --report "{output.report}" \
            --keep-files "{params.keep_files}" \
            --remove-files "{params.remove_files}"
        """

rule apply_external_sample_set:
    input:
        bed=rules.standardize_genotype.output.bed,
        bim=rules.standardize_genotype.output.bim,
        fam=rules.standardize_genotype.output.fam,
        keep=rules.build_external_sample_set.output.keep,
    output:
        bed=ws_path(
            "external_samples/genotype.bed"
        ),
        bim=ws_path(
            "external_samples/genotype.bim"
        ),
        fam=ws_path(
            "external_samples/genotype.fam"
        ),
        log=ws_path(
            "external_samples/genotype.log"
        ),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=90,
        mem_mb=12000,
    params:
        source=ws_path("standardization/genotype"),
        prefix=ws_path("external_samples/genotype"),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.bed}")"

        plink2 \
            --bfile "{params.source}" \
            --keep "{input.keep}" \
            --make-bed \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory {resources.mem_mb}
        """

SAMPLE_MISSINGNESS_PREFIX = ("sample_missingness/sample_missingness")

rule sample_missingness_report:
    input:
        # Standardized genotype dataset after applying external
        # sample inclusion/exclusion lists.
        bed=rules.apply_external_sample_set.output.bed,
        bim=rules.apply_external_sample_set.output.bim,
        fam=rules.apply_external_sample_set.output.fam,
    output:
        # Per-sample genotype missingness statistics.
        smiss=ws_path(
            "sample_missingness/sample_missingness.smiss"
        ),

        # PLINK execution log.
        log=ws_path(
            "sample_missingness/sample_missingness.log"
        ),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=60,
        mem_mb=8000,
    params:
        # Prefix of the PLINK dataset produced by
        # apply_external_sample_set.
        bfile=ws_path(
            "external_samples/genotype"
        ),

        # Prefix used for the PLINK output files.
        prefix=ws_path(
            "sample_missingness/sample_missingness"
        ),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.smiss}")"

        plink2 \
            --bfile "{params.bfile}" \
            --autosome \
            --missing sample-only \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory {resources.mem_mb}
        """

rule plot_sample_missingness:
    input:
        smiss=rules.sample_missingness_report.output.smiss,
    output:
        png=ws_path(
            "sample_missingness/sample_missingness.png"
        ),
        summary=ws_path(
            "sample_missingness/sample_missingness.summary.tsv"
        ),
    conda:
        "../envs/r_environment.yaml"
    params:
        threshold=lambda wc: config.get(
            "sample_qc", {}
        ).get("mind", 0.1),
    shell:
        r"""
        Rscript workflow/scripts/plot_sample_missingness.R \
            --smiss {input.smiss} \
            --threshold {params.threshold} \
            --plot {output.png} \
            --summary {output.summary}
        """
