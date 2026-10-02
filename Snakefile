configfile: "config/config.yaml"


# Load shared helper functions.
include: "workflow/rules/common.smk"


# Sex-QC -----------------------------------------------------------------------

# Diagnostic outputs for reviewing sex-QC thresholds.
SEX_REVIEW_OUTPUTS = [
    ws_path("review/sex/sex_F_distribution.pdf"),
    ws_path("review/sex/sex_candidate_classification.tsv"),
    ws_path("review/sex/candidate_sex_thresholds.yaml"),
    ws_path("review/sex/sex_qc.summary.tsv"),
]

# Final sex-QC outputs generated using the currently configured thresholds.
SEX_FINAL_OUTPUTS = [
    ws_path("sex_qc/automatic_sex_exclusions.tsv"),
    ws_path("sex_qc/sex_classification.tsv"),
    ws_path("sex_qc/accepted_sex_qc.summary.tsv"),
]

# Always generate both the diagnostic review outputs and the final
# sex-QC classification/exclusion outputs.
SEX_QC_OUTPUTS = SEX_REVIEW_OUTPUTS + SEX_FINAL_OUTPUTS


# Population-structure and relatedness QC -------------------------------------

# These outputs are always part of the final workflow.
# The structure/KING branch therefore runs during the first execution using
# the currently configured sex thresholds. If the thresholds are revised
# later, Snakemake will rerun the sex-dependent downstream rules.
STRUCTURE_OUTPUTS = [
    ws_path("structure/structure_samples.keep"),
    ws_path("structure/structure_sample_exclusions.tsv"),
    ws_path("structure/structure_samples.summary.tsv"),

    ws_path("structure/genotype.bed"),
    ws_path("structure/genotype.bim"),
    ws_path("structure/genotype.fam"),
    ws_path("structure/genotype.log"),

    ws_path("structure/pcadapt/king_variants.keep.txt"),
    ws_path("structure/pcadapt/ancestry_associated_variants.tsv"),
    ws_path("structure/pcadapt/pcadapt_diagnostics.pdf"),
    ws_path("structure/pcadapt/pcadapt.summary.tsv"),

    ws_path("structure/king.king.cutoff.in.id"),
    ws_path("structure/king.king.cutoff.out.id"),
    ws_path("structure/king.log"),

    ws_path("structure/king_pairs.kin0"),
    ws_path("structure/king_pairs.log"),
]


rule all:
    input:
        [
            # Input inspection
            ws_path("standardization/input.ok"),
            ws_path("standardization/input.summary.tsv"),

            # Select genotype samples before matching reported sex
            ws_path("external_samples/eligible_samples.keep"),
            ws_path("external_samples/external_sample_exclusions.tsv"),
            ws_path("external_samples/external_samples.report.tsv"),
            ws_path("external_samples/genotype.bed"),
            ws_path("external_samples/genotype.bim"),
            ws_path("external_samples/genotype.fam"),

            # Reported sex
            ws_path("standardization/reported_sex.update.tsv"),
            ws_path("standardization/phenotype_matching_report.tsv"),
            ws_path("standardization/genotype_without_phenotype.tsv"),
            ws_path("standardization/phenotype_without_genotype.tsv"),

            # Standardized genotype dataset
            ws_path("standardization/genotype.bed"),
            ws_path("standardization/genotype.bim"),
            ws_path("standardization/genotype.fam"),

            # Sample missingness
            ws_path("sample_missingness/sample_missingness.smiss"),
            ws_path("sample_missingness/sample_missingness.log"),
            ws_path("sample_missingness/sample_missingness.png"),
            ws_path("sample_missingness/sample_missingness.summary.tsv"),
            ws_path("sample_missingness/passing_samples.id"),
            ws_path("sample_missingness/mind_filter.log"),

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

            # LD-pruned X markers used for the candidate sex check
            ws_path("sex_qc/x_pruned/markers.prune.in"),
            ws_path("sex_qc/x_pruned/markers.prune.out"),
            ws_path("sex_qc/x_pruned/markers.log"),

            # Candidate sex check
            ws_path("sex_qc/sex_candidate.sexcheck"),
            ws_path("sex_qc/sex_candidate.log"),

            # Sex-QC diagnostic and final outputs
            *SEX_QC_OUTPUTS,

            # Population structure and KING relatedness QC
            *STRUCTURE_OUTPUTS,
        ]


# Include workflow modules.
include: "workflow/rules/standardize_input_data.smk"
include: "workflow/rules/sample_missingness.smk"
include: "workflow/rules/autosomal_qc.smk"
include: "workflow/rules/heterozygosity_qc.smk"
include: "workflow/rules/sex_qc.smk"
include: "workflow/rules/structure.smk"
