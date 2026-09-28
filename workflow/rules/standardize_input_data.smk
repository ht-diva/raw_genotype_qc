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
        # Work only with samples selected at the start of the pipeline.
        fam=ws_path("external_samples/genotype.fam"),

        # Phenotype containing sample IDs and reported sex.
        phenotype=config["phenotype"]["path"],

        # Ensures that input validation is completed successfully
        # before this rule starts.
        validated=rules.inspect_bed_input.output.ok,
    output:
        # PLINK-compatible file containing: FID, IID, SEX.
        update_sex=ws_path(
            "standardization/reported_sex.update.tsv"
        ),

        # Summary of genotype-phenotype matching.
        report=ws_path(
            "standardization/phenotype_matching_report.tsv"
        ),

        # Genotyped samples that were not found in the phenotype.
        unmatched_genotypes=ws_path(
            "standardization/genotype_without_phenotype.tsv"
        ),

        # Phenotype records that were not found in the genotype data.
        unmatched_phenotype=ws_path(
            "standardization/phenotype_without_genotype.tsv"
        ),

    params:
        # Phenotype column containing the sample identifier.
        phenotype_id_col=lambda wc: cfg(
            "phenotype/phenotype_id_col",
            "IID",
        ),

        # Phenotype column containing reported sex.
        sex_col=lambda wc: cfg(
            "phenotype/sex_col",
            "SEX",
        ),

        # Method used to match FAM and phenotype sample IDs.
        genotype_id_mode=lambda wc: cfg(
            "phenotype/genotype_id_mode",
            "direct",
        ),

        # Regular expression used only in regex matching mode.
        genotype_id_regex=lambda wc: (
            cfg("phenotype/genotype_id_regex") or "NA"
        ),

        # Method used to normalize sample IDs before matching.
        id_normalization=lambda wc: cfg(
            "phenotype/id_normalization",
            "string",
        ),

        # Phenotype values interpreted as male.
        male_values=lambda wc: ",".join(
            map(
                str,
                cfg(
                    "phenotype/male_values",
                    [1, "M", "Male"],
                ),
            )
        ),

        # Phenotype values interpreted as female.
        female_values=lambda wc: ",".join(
            map(
                str,
                cfg(
                    "phenotype/female_values",
                    [2, "F", "Female"],
                ),
            )
        ),
    conda:
        "../envs/r_environment.yaml"
    shell:
        r"""
        set -euo pipefail

        Rscript workflow/scripts/reported_sex.R \
            --fam "{input.fam}" \
            --phenotype "{input.phenotype}" \
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
            --unmatched-phenotype "{output.unmatched_phenotype}"
        """


rule standardize_genotype_pgen:
    input:
        bed=ws_path("external_samples/genotype.bed"),
        bim=ws_path("external_samples/genotype.bim"),
        fam=ws_path("external_samples/genotype.fam"),
        chromosome_map=rules.inspect_bed_input.output.chr_map,
        update_sex=rules.prepare_reported_sex.output.update_sex,
    output:
        pgen=temp(
            ws_path("standardization/genotype_sorted.pgen")
        ),
        pvar=temp(
            ws_path("standardization/genotype_sorted.pvar")
        ),
        psam=temp(
            ws_path("standardization/genotype_sorted.psam")
        ),
        log=temp(
            ws_path("standardization/genotype_sorted.log")
        ),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=90,
        mem_mb=12000,
    params:
        source=ws_path("external_samples/genotype"),
        prefix=ws_path(
            "standardization/genotype_sorted"
        ),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.pgen}")"

        plink2 \
            --bfile "{params.source}" \
            --rename-chrs "{input.chromosome_map}" \
            --sort-vars natural \
            --update-sex "{input.update_sex}" \
            --make-pgen \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory {resources.mem_mb}
        """


rule standardize_genotype:
    input:
        pgen=rules.standardize_genotype_pgen.output.pgen,
        pvar=rules.standardize_genotype_pgen.output.pvar,
        psam=rules.standardize_genotype_pgen.output.psam,
    output:
        bed=ws_path(
            "standardization/genotype.bed"
        ),
        bim=ws_path(
            "standardization/genotype.bim"
        ),
        fam=ws_path(
            "standardization/genotype.fam"
        ),
        log=ws_path(
            "standardization/genotype.log"
        ),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=90,
        mem_mb=12000,
    params:
        source=ws_path(
            "standardization/genotype_sorted"
        ),
        prefix=ws_path(
            "standardization/genotype"
        ),
    shell:
        r"""
        set -euo pipefail

        plink2 \
            --pfile "{params.source}" \
            --make-bed \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory {resources.mem_mb}
        """
