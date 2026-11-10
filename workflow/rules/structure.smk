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
        1
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
            --make-king-table cols=fid,id,nsnp,hethet,ibs0,kinship \
            --king-table-filter {params.cutoff} \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory {resources.mem_mb}
        """

rule king_unrelated:
    input:
        bed=rules.create_structure_qc_set.output.bed,
        bim=rules.create_structure_qc_set.output.bim,
        fam=rules.create_structure_qc_set.output.fam,
        kin0=rules.king_table.output.kin0,
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
            --king-cutoff-table "{input.kin0}" {params.cutoff} \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory {resources.mem_mb}
        """




# Reference inputs are required only when projection is explicitly enabled.
def external_reference_inputs(wildcards):
    prefix = cfg("pca_reference/bed_prefix", None)
    labels = cfg("pca_reference/labels_path", None)
    if not prefix or not labels:
        raise WorkflowError(
            "Reference projection requires pca_reference/bed_prefix and labels_path."
        )
    return [str(prefix) + suffix for suffix in (".bed", ".bim", ".fam")] + [labels]


rule pca_initial:
    input:
        script="workflow/scripts/pca_diagnostics.R",
        helpers="workflow/scripts/pca_helpers.R",
        bed=rules.create_structure_qc_set.output.bed,
        bim=rules.create_structure_qc_set.output.bim,
        fam=rules.create_structure_qc_set.output.fam,
        unrelated=rules.king_unrelated.output.keep,
        keep=rules.structure_sample_set.output.keep,
    output:
        pcs=ws_path("review/pca/pca_initial_all_samples.tsv"),
        scores_plot=ws_path("review/pca/pca_scores.pdf"),
        loadings_plot=ws_path("review/pca/pca_loadings.pdf"),
        loadings=ws_path("review/pca/pca_loadings.tsv"),
        lrldr=ws_path("review/pca/autoSVD_long_range_LD_regions.tsv"),
        robust=ws_path("review/pca/robust_pca_distances.tsv"),
        suggested=ws_path("review/pca/suggested_pca_outliers.tsv"),
        retained=ws_path("review/pca/autoSVD_retained_variants.txt"),
        summary=ws_path("review/pca/pca_initial_summary.tsv"),
    conda:
        "../envs/r_environment.yaml"
    threads: 12
    resources:
        runtime=360,
        mem_mb=65536,
    params:
        prefix=ws_path(STRUCTURE_QC_PREFIX),
        # Accept the supplied rule's key and the existing pipeline's key.
        n_pcs=lambda wc: cfg("pca/n_pcs", cfg("pca/components", 20)),
        robust_pcs=lambda wc: cfg("pca/robust_distance_n_pcs", 6),
        robust_mad=lambda wc: cfg("pca/robust_distance_mad_multiplier", 6),
        max_iter=lambda wc: cfg("pca/autosvd_max_iter", 5),
        min_maf=lambda wc: cfg("pca/autosvd_min_maf", 0.01),
        seed=lambda wc: cfg("pca/seed", 1),
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.pcs:q})"
        Rscript {input.script:q} \
          --bed-prefix {params.prefix:q} --unrelated-file {input.unrelated:q} \
          --keep-file {input.keep:q} --n-pcs {params.n_pcs} \
          --robust-distance-n-pcs {params.robust_pcs} \
          --robust-mad-multiplier {params.robust_mad} \
          --autosvd-max-iter {params.max_iter} \
          --autosvd-min-maf {params.min_maf} --seed {params.seed} \
          --pcs-output {output.pcs:q} --scores-plot {output.scores_plot:q} \
          --loadings-plot {output.loadings_plot:q} --loadings-output {output.loadings:q} \
          --lrldr-output {output.lrldr:q} --robust-output {output.robust:q} \
          --suggested-outliers {output.suggested:q} \
          --retained-variants {output.retained:q} \
          --summary-output {output.summary:q} --threads {threads}
        '''


if cfg("pca_reference/enabled", False):
    rule external_reference_projection:
        input:
            script="workflow/scripts/external_reference_projection.R",
            helpers="workflow/scripts/pca_helpers.R",
            study_bed=rules.create_structure_qc_set.output.bed,
            study_bim=rules.create_structure_qc_set.output.bim,
            study_fam=rules.create_structure_qc_set.output.fam,
            study_pcs=rules.pca_initial.output.pcs,
            reference=external_reference_inputs,
        output:
            assignments=ws_path("review/pca/external_reference_assignments.tsv"),
            plot=ws_path("review/pca/study_pca_colored_by_reference.pdf"),
            summary=ws_path("review/pca/external_reference.summary.tsv"),
        conda:
            "../envs/r_environment.yaml"
        threads: 8
        resources:
            runtime=360,
            mem_mb=65536,
        params:
            enabled="True",
            study_prefix=ws_path(STRUCTURE_QC_PREFIX),
            ref_prefix=lambda wc: cfg("pca_reference/bed_prefix", None) or "NA",
            labels=lambda wc: cfg("pca_reference/labels_path", None) or "NA",
            seed=lambda wc: cfg("pca_reference/seed", 1),
            fid_col=lambda wc: cfg("pca_reference/fid_col", "FID"),
            iid_col=lambda wc: cfg("pca_reference/iid_col", "IID"),
            group_col=lambda wc: cfg("pca_reference/group_col", "GROUP"),
            build_study=lambda wc: cfg("pca_reference/study_build", cfg("genome_build", "b38")),
            build_ref=lambda wc: cfg("pca_reference/reference_build", cfg("genome_build", "b38")),
            n_pcs=lambda wc: cfg("pca_reference/n_pcs", 10),
            match_min=lambda wc: cfg("pca_reference/match_min_prop", 0.5),
            quantile=lambda wc: cfg("pca_reference/assignment_quantile", 0.995),
            multiplier=lambda wc: cfg("pca_reference/assignment_radius_multiplier", 1.2),
        shell:
            r'''
            set -euo pipefail
            mkdir -p "$(dirname {output.assignments:q})"
            Rscript {input.script:q} \
              --enabled {params.enabled} --study-prefix {params.study_prefix:q} \
              --study-pcs {input.study_pcs:q} --reference-prefix {params.ref_prefix:q} \
              --labels {params.labels:q} --fid-col {params.fid_col:q} --iid-col {params.iid_col:q} \
              --group-col {params.group_col:q} --study-build {params.build_study:q} \
              --reference-build {params.build_ref:q} --n-pcs {params.n_pcs} \
              --match-min-prop {params.match_min} \
              --assignment-quantile {params.quantile} \
              --assignment-radius-multiplier {params.multiplier} --seed {params.seed} \
              --assignments {output.assignments:q} --plot-output {output.plot:q} \
              --summary-output {output.summary:q} --threads {threads}
            '''




rule pca_review_bundle:
    input:
        diagnostics=rules.pca_initial.output,
        kinship=rules.king_table.output.kin0,
        reference=lambda wc: (
            [ws_path("review/pca/external_reference_assignments.tsv"),
             ws_path("review/pca/study_pca_colored_by_reference.pdf"),
             ws_path("review/pca/external_reference.summary.tsv")]
            if cfg("pca_reference/enabled", False) else []
        ),
    output:
        ready=touch(ws_path("review/pca/REVIEW_REQUIRED.done")),
    shell:
        r"""
        set -euo pipefail
        echo "PCA review required: inspect scores, loadings, long-range LD regions, robust distances, suggested outliers and optional reference assignments. This marker records completed diagnostics, not accepted exclusions. No suggested PCA outliers are removed automatically." >&2
        touch {output.ready:q}
        """
