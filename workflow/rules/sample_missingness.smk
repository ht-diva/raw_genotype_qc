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

