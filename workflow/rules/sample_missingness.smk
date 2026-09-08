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


