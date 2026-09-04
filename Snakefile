configfile: "config/config.yaml"


# Shared helpers must be loaded before any rule file.
include: "workflow/rules/common.smk"


rule all:
    input:
        ws_path("standardization/input.ok"),
        ws_path("standardization/reported_sex.update.tsv"),
        ws_path("standardization/sample_metadata_report.tsv"),
        ws_path("standardization/genotype_without_metadata.tsv"),
        ws_path("standardization/metadata_without_genotype.tsv"),


# Include only modules which are currently implemented.
include: "workflow/rules/standardize_input_data.smk"
