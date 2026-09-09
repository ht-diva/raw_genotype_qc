configfile: "config/config.yaml"


# Load shared helper functions.
include: "workflow/rules/common.smk"


rule all:
    input:
        [
            # Input inspection
            ws_path("standardization/input.ok"),
            ws_path("standardization/input.summary.tsv"),

            # Reported sex
            ws_path("standardization/reported_sex.update.tsv"),
            ws_path("standardization/sample_metadata_report.tsv"),
            ws_path("standardization/genotype_without_metadata.tsv"),
            ws_path("standardization/metadata_without_genotype.tsv"),

            # Standardized genotype dataset
            ws_path("standardization/genotype.bed"),
            ws_path("standardization/genotype.bim"),
            ws_path("standardization/genotype.fam"),

            # External sample selection
            ws_path("external_samples/eligible_samples.keep"),
            ws_path(
                "external_samples/external_sample_exclusions.tsv"
            ),
            ws_path(
                "external_samples/external_samples.report.tsv"
            ),

            # Dataset after external sample selection
            ws_path("external_samples/genotype.bed"),
            ws_path("external_samples/genotype.bim"),
            ws_path("external_samples/genotype.fam"),

            # Sample missingness
            ws_path(
                "sample_missingness/sample_missingness.smiss"
            ),
            ws_path(
                "sample_missingness/sample_missingness.log"
            ),
            ws_path(
                "sample_missingness/sample_missingness.png"
            ),
            ws_path(
                "sample_missingness/sample_missingness.summary.tsv"
            ),
            ws_path(
                "sample_missingness/passing_samples.id"
            ),
            ws_path(
                "sample_missingness/mind_filter.log"
            ),

            # Final genotype dataset after sample missingness filtering
            ws_path("sample_missingness/genotype.bed"),
            ws_path("sample_missingness/genotype.bim"),
            ws_path("sample_missingness/genotype.fam"),
            ws_path("sample_missingness/genotype.log"),
        ]


# Include workflow modules.
include: "workflow/rules/standardize_input_data.smk"
include: "workflow/rules/sample_missingness.smk"