"""Population-structure filtering and relatedness QC."""


STRUCTURE_QC_PREFIX = "structure/genotype"


rule structure_sample_set:
    input:
        fam=rules.apply_heterozygosity_exclusions.output.fam,
        sex_exclusions=rules.apply_sex_thresholds.output.exclusions,
    output:
        keep=ws_path("structure/structure_samples.keep"),
        exclusions=ws_path("structure/structure_sample_exclusions.tsv"),
        summary=ws_path("structure/structure_samples.summary.tsv"),
    conda:
        "../envs/r_environment.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.keep}")"

        Rscript workflow/scripts/build_sample_set.R \
            --fam "{input.fam}" \
            --exclude-files "{input.sex_exclusions}" \
            --exclusions-output "{output.exclusions}" \
            --keep-output "{output.keep}" \
            --summary "{output.summary}"
        """


rule create_structure_qc_set:
    input:
        bed=rules.apply_heterozygosity_exclusions.output.bed,
        bim=rules.apply_heterozygosity_exclusions.output.bim,
        fam=rules.apply_heterozygosity_exclusions.output.fam,
        keep=rules.structure_sample_set.output.keep,
    output:
        bed=ws_path(STRUCTURE_QC_PREFIX + ".bed"),
        bim=ws_path(STRUCTURE_QC_PREFIX + ".bim"),
        fam=ws_path(STRUCTURE_QC_PREFIX + ".fam"),
        log=ws_path(STRUCTURE_QC_PREFIX + ".log"),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=120,
        mem_mb=16384,
    params:
        bfile=ws_path("heterozygosity/genotype"),
        prefix=ws_path(STRUCTURE_QC_PREFIX),
        maf=lambda wc: cfg("thresholds/maf"),
        geno=lambda wc: cfg("thresholds/variant_missingness"),
        hwe=lambda wc: cfg("thresholds/hwe"),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.bed}")"

        plink2 \
            --bfile "{params.bfile}" \
            --keep "{input.keep}" \
            --autosome \
            --maf "{params.maf}" \
            --geno "{params.geno}" \
            --hwe "{params.hwe}" 0 \
            --make-bed \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory {resources.mem_mb}
        """


rule pcadapt_for_king:
    input:
        bed=rules.create_structure_qc_set.output.bed,
        bim=rules.create_structure_qc_set.output.bim,
        fam=rules.create_structure_qc_set.output.fam,
    output:
        keep_variants=ws_path("structure/pcadapt/king_variants.keep.txt"),
        excluded_variants=ws_path(
            "structure/pcadapt/ancestry_associated_variants.tsv"
        ),
        plot=ws_path("structure/pcadapt/pcadapt_diagnostics.pdf"),
        summary=ws_path("structure/pcadapt/pcadapt.summary.tsv"),
    conda:
        "../envs/r_environment.yaml"
    threads:
        8
    resources:
        runtime=240,
        mem_mb=65536,
    params:
        prefix=ws_path(STRUCTURE_QC_PREFIX),
        n_pcs=lambda wc: cfg("pcadapt/n_pcs", 5),
        max_iter=lambda wc: cfg(
            "pcadapt/preliminary_autosvd_max_iter",
            1,
        ),
        min_maf=lambda wc: cfg("pcadapt/autosvd_min_maf", 0.01),
        keep_p=lambda wc: cfg(
            "pcadapt/keep_adjusted_p_above",
            0.05,
        ),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.keep_variants}")"

        Rscript workflow/scripts/pcadapt_for_king.R \
            --bed-prefix "{params.prefix}" \
            --n-pcs {params.n_pcs} \
            --autosvd-max-iter {params.max_iter} \
            --autosvd-min-maf {params.min_maf} \
            --keep-adjusted-p-above {params.keep_p} \
            --keep-variants "{output.keep_variants}" \
            --excluded-variants "{output.excluded_variants}" \
            --plot-output "{output.plot}" \
            --summary-output "{output.summary}" \
            --threads {threads}
        """


rule king_unrelated:
    input:
        bed=rules.create_structure_qc_set.output.bed,
        bim=rules.create_structure_qc_set.output.bim,
        fam=rules.create_structure_qc_set.output.fam,
        variants=rules.pcadapt_for_king.output.keep_variants,
    output:
        keep=ws_path("structure/king.king.cutoff.in.id"),
        removed=ws_path("structure/king.king.cutoff.out.id"),
        log=ws_path("structure/king.log"),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        12
    resources:
        runtime=360,
        mem_mb=49152,
    params:
        bfile=ws_path(STRUCTURE_QC_PREFIX),
        prefix=ws_path("structure/king"),
        cutoff=lambda wc: cfg("thresholds/king_cutoff"),
    shell:
        r"""
        set -euo pipefail

        plink2 \
            --bfile "{params.bfile}" \
            --extract "{input.variants}" \
            --king-cutoff {params.cutoff} \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory {resources.mem_mb}
        """


rule king_table:
    input:
        bed=rules.create_structure_qc_set.output.bed,
        bim=rules.create_structure_qc_set.output.bim,
        fam=rules.create_structure_qc_set.output.fam,
        variants=rules.pcadapt_for_king.output.keep_variants,
    output:
        kin0=ws_path("structure/king_pairs.kin0"),
        log=ws_path("structure/king_pairs.log"),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        12
    resources:
        runtime=360,
        mem_mb=49152,
    params:
        bfile=ws_path(STRUCTURE_QC_PREFIX),
        prefix=ws_path("structure/king_pairs"),
        cutoff=lambda wc: cfg("thresholds/king_cutoff"),
    shell:
        r"""
        set -euo pipefail

        plink2 \
            --bfile "{params.bfile}" \
            --extract "{input.variants}" \
            --make-king-table \
            --king-table-filter {params.cutoff} \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory {resources.mem_mb}
        """