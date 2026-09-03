rule inspect_merged_input:
    input:
        bed=rules.merge_canonical_chromosomes.output.bed,
        bim=rules.merge_canonical_chromosomes.output.bim,
        fam=rules.merge_canonical_chromosomes.output.fam,
    output:
        chr_map=ws_path("01_standardization/chromosome_map.tsv"),
        summary=ws_path("01_standardization/input.summary.tsv"),
        duplicate_samples=ws_path("01_standardization/duplicate_samples.tsv"),
        duplicate_variant_ids=ws_path("01_standardization/duplicate_variant_ids.tsv"),
        duplicate_positions=ws_path("01_standardization/duplicate_positions.tsv"),
        ok=touch(ws_path("01_standardization/input.ok")),
    params:
        fail_duplicate_samples=lambda wc: str(cfg("input_validation/fail_on_duplicate_samples", True)),
    shell:
        r'''
        mkdir -p $(dirname {output.summary})
        workflow/scripts/inspect_bed_input.sh \
          --bed {input.bed} --bim {input.bim} --fam {input.fam} \
          --chromosome-map {output.chr_map} \
          --summary {output.summary} \
          --duplicate-samples {output.duplicate_samples} \
          --duplicate-variant-ids {output.duplicate_variant_ids} \
          --duplicate-positions {output.duplicate_positions} \
          --fail-duplicate-samples {params.fail_duplicate_samples}
        touch {output.ok}
        '''


rule prepare_reported_sex:
    input:
        fam=rules.merge_canonical_chromosomes.output.fam,
        validated=rules.inspect_merged_input.output.ok,
    output:
        update_sex=ws_path("01_standardization/reported_sex.update.tsv"),
        report=ws_path("01_standardization/sample_metadata_report.tsv"),
        unmatched_genotypes=ws_path("01_standardization/genotype_without_metadata.tsv"),
        unmatched_metadata=ws_path("01_standardization/metadata_without_genotype.tsv"),
    conda:
        "envs/r_qc.yaml"
    params:
        metadata_path=lambda wc: cfg("sample_metadata/path") or "NA",
        phenotype_id_col=lambda wc: cfg("sample_metadata/phenotype_id_col", "IID"),
        sex_col=lambda wc: cfg("sample_metadata/sex_col", "SEX"),
        genotype_id_mode=lambda wc: cfg("sample_metadata/genotype_id_mode", "fam"),
        genotype_id_regex=lambda wc: cfg("sample_metadata/genotype_id_regex") or "NA",
        id_normalization=lambda wc: cfg("sample_metadata/id_normalization", "string"),
        male_values=lambda wc: ",".join(map(str, cfg("sample_metadata/male_values", [1, "M", "Male"]))),
        female_values=lambda wc: ",".join(map(str, cfg("sample_metadata/female_values", [2, "F", "Female"]))),
    shell:
        r'''
        Rscript workflow/scripts/prepare_reported_sex.R \
          --fam {input.fam} \
          --metadata '{params.metadata_path}' \
          --phenotype-id-col '{params.phenotype_id_col}' \
          --sex-col '{params.sex_col}' \
          --genotype-id-mode {params.genotype_id_mode} \
          --genotype-id-regex '{params.genotype_id_regex}' \
          --id-normalization {params.id_normalization} \
          --male-values '{params.male_values}' \
          --female-values '{params.female_values}' \
          --update-sex {output.update_sex} \
          --report {output.report} \
          --unmatched-genotypes {output.unmatched_genotypes} \
          --unmatched-metadata {output.unmatched_metadata}
        '''


rule standardize_genotype:
    input:
        bed=rules.merge_canonical_chromosomes.output.bed,
        bim=rules.merge_canonical_chromosomes.output.bim,
        fam=rules.merge_canonical_chromosomes.output.fam,
        chr_map=rules.inspect_merged_input.output.chr_map,
        update_sex=rules.prepare_reported_sex.output.update_sex,
    output:
        bed=ws_path(STANDARD_PREFIX + ".bed"),
        bim=ws_path(STANDARD_PREFIX + ".bim"),
        fam=ws_path(STANDARD_PREFIX + ".fam"),
        log=ws_path(STANDARD_PREFIX + ".log"),
    threads: 8
    resources:
        runtime=120,
        mem_mb=16384,
    params:
        plink2=PLINK2,
        source=ws_path(MERGED_INPUT_PREFIX),
        prefix=ws_path(STANDARD_PREFIX),
    shell:
        r'''
        {params.plink2} \
          --bfile {params.source} \
          --rename-chrs {input.chr_map} \
          --update-sex {input.update_sex} \
          --sort-vars \
          --make-bed \
          --out {params.prefix} \
          --threads {threads} --memory 16000
        '''