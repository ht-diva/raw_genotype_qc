from pathlib import Path


configfile: "config/config.yaml"


# Load shared helper functions.
include: "workflow/rules/common.smk"


# Sex-QC manual approval --------------------------------------------------------

SEX_APPROVAL_FILE = cfg(
    "accepted_sex_threshold",
    "config/accepted_sex_threshold.yaml",
)

SEX_REVIEW_OUTPUTS = [
    ws_path("review/sex/sex_F_distribution.pdf"),
    ws_path("review/sex/sex_candidate_classification.tsv"),
    ws_path("review/sex/candidate_sex_threshold.yaml"),
    ws_path("review/sex/sex_qc.summary.tsv"),
]

SEX_FINAL_OUTPUTS = [
    ws_path("sex_qc/automatic_sex_exclusions.tsv"),
    ws_path("sex_qc/sex_classification.tsv"),
    ws_path("sex_qc/accepted_sex_qc.summary.tsv"),
]

SEX_QC_OUTPUTS = (
    SEX_FINAL_OUTPUTS
    if Path(SEX_APPROVAL_FILE).is_file()
    else SEX_REVIEW_OUTPUTS
)


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

            # Autosomal dataset for sample QC
            ws_path("autosomal_qc/genotype.bed"),
            ws_path("autosomal_qc/genotype.bim"),
            ws_path("autosomal_qc/genotype.fam"),
            ws_path("autosomal_qc/genotype.log"),

            # Heterozygosity QC
            ws_path("heterozygosity/heterozygosity.het"),
            ws_path("heterozygosity/heterozygosity.log"),
            ws_path("heterozygosity/heterozygosity.png"),
            ws_path("heterozygosity/heterozygosity_exclusions.tsv"),
            ws_path("heterozygosity/heterozygosity_all_samples.tsv"),
            ws_path("heterozygosity/heterozygosity.summary.tsv"),

            # Dataset after heterozygosity filtering
            ws_path("heterozygosity/genotype.bed"),
            ws_path("heterozygosity/genotype.bim"),
            ws_path("heterozygosity/genotype.fam"),
            ws_path("heterozygosity/genotype.log"),

            # Sex QC dataset after PAR splitting
            ws_path("sex_qc/split_par/genotype.bed"),
            ws_path("sex_qc/split_par/genotype.bim"),
            ws_path("sex_qc/split_par/genotype.fam"),
            ws_path("sex_qc/split_par/genotype.log"),

            # X-chromosome markers used for sex QC
            ws_path("sex_qc/x_qc/genotype.bed"),
            ws_path("sex_qc/x_qc/genotype.bim"),
            ws_path("sex_qc/x_qc/genotype.fam"),
            ws_path("sex_qc/x_qc/genotype.log"),

            # Candidate sex check
            ws_path("sex_qc/sex_candidate.sexcheck"),
            ws_path("sex_qc/sex_candidate.log"),

            # Sex-QC review or final outputs
            *SEX_QC_OUTPUTS,
        ]


# Include workflow modules.
include: "workflow/rules/standardize_input_data.smk"
include: "workflow/rules/sample_missingness.smk"
include: "workflow/rules/autosomal_qc.smk"
include: "workflow/rules/heterozygosity_qc.smk"
include: "workflow/rules/sex_qc.smk"