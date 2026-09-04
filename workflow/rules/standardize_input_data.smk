rule inspect_bed_input:
    input:
        bed=config["bed_path"],
        bim=config["bim_path"],
        fam=config["fam_path"],
    output:
        chr_map=ws_path(
            "standardization/chromosome_map.tsv"
        ),
        summary=ws_path(
            "standardization/input.summary.tsv"
        ),
        duplicate_samples=ws_path(
            "standardization/duplicate_samples.tsv"
        ),
        duplicate_variant_ids=ws_path(
            "standardization/duplicate_variant_ids.tsv"
        ),
        duplicate_positions=ws_path(
            "standardization/duplicate_positions.tsv"
        ),
        ok=touch(
            ws_path("standardization/input.ok")
        ),
    params:
        fail_duplicate_samples=lambda wc: str(
            cfg(
                "input_validation/fail_on_duplicate_samples",
                True,
            )
        ).lower(),
    shell:
        r"""
        set -euo pipefail

        bash workflow/scripts/inspect_bed_input.sh \
            --bed "{input.bed}" \
            --bim "{input.bim}" \
            --fam "{input.fam}" \
            --chromosome-map "{output.chr_map}" \
            --summary "{output.summary}" \
            --duplicate-samples "{output.duplicate_samples}" \
            --duplicate-variant-ids "{output.duplicate_variant_ids}" \
            --duplicate-positions "{output.duplicate_positions}" \
            --fail-duplicate-samples \
                "{params.fail_duplicate_samples}"
        """

rule prepare_reported_sex:
    input:
        # FAM file already validated by inspect_bed_input.
        fam=config["fam_path"],

        # Metadata containing sample IDs and reported sex.
        metadata=config["sample_metadata"]["path"],

        # Ensures that input validation is completed successfully
        # before this rule starts.
        validated=rules.inspect_bed_input.output.ok,
    output:
        # PLINK-compatible file containing: FID, IID, SEX.
        update_sex=ws_path(
            "standardization/reported_sex.update.tsv"
        ),

        # Summary of genotype-metadata matching.
        report=ws_path(
            "standardization/sample_metadata_report.tsv"
        ),

        # Genotyped samples that were not found in the metadata.
        unmatched_genotypes=ws_path(
            "standardization/genotype_without_metadata.tsv"
        ),

        # Metadata records that were not found in the genotype data.
        unmatched_metadata=ws_path(
            "standardization/metadata_without_genotype.tsv"
        ),

    params:
        # Metadata column containing the sample identifier.
        phenotype_id_col=lambda wc: cfg(
            "sample_metadata/phenotype_id_col",
            "IID",
        ),

        # Metadata column containing reported sex.
        sex_col=lambda wc: cfg(
            "sample_metadata/sex_col",
            "SEX",
        ),

        # Method used to match FAM and metadata sample IDs.
        genotype_id_mode=lambda wc: cfg(
            "sample_metadata/genotype_id_mode",
            "direct",
        ),

        # Regular expression used only in regex matching mode.
        genotype_id_regex=lambda wc: (
            cfg("sample_metadata/genotype_id_regex") or "NA"
        ),

        # Method used to normalize sample IDs before matching.
        id_normalization=lambda wc: cfg(
            "sample_metadata/id_normalization",
            "string",
        ),

        # Metadata values interpreted as male.
        male_values=lambda wc: ",".join(
            map(
                str,
                cfg(
                    "sample_metadata/male_values",
                    [1, "M", "Male"],
                ),
            )
        ),

        # Metadata values interpreted as female.
        female_values=lambda wc: ",".join(
            map(
                str,
                cfg(
                    "sample_metadata/female_values",
                    [2, "F", "Female"],
                ),
            )
        ),
    shell:
        r"""
        set -euo pipefail

        Rscript workflow/scripts/prepare_reported_sex.R \
            --fam "{input.fam}" \
            --metadata "{input.metadata}" \
            --phenotype-id-col "{params.phenotype_id_col}" \
            --sex-col "{params.sex_col}" \
            --genotype-id-mode "{params.genotype_id_mode}" \
            --genotype-id-regex "{params.genotype_id_regex}" \
            --id-normalization "{params.id_normalization}" \
            --male-values "{params.male_values}" \
            --female-values "{params.female_values}" \
            --update-sex "{output.update_sex}" \
            --report "{output.report}" \
            --unmatched-genotypes "{output.unmatched_genotypes}" \
            --unmatched-metadata "{output.unmatched_metadata}"
        """