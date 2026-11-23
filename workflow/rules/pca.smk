"""Study PCA after KING and optional external-reference review."""


PCA_METHOD = cfg("pca/method", "autosvd")
if PCA_METHOD not in ("plink2", "autosvd"):
    raise WorkflowError("pca/method must be 'plink2' or 'autosvd'.")

PCA_FIT_PREFIX = ws_path("structure/pca/unrelated")
PCA_SCORE_PREFIX = ws_path("structure/pca/all_samples")
PCA_COMPONENTS = cfg("pca/n_pcs", cfg("pca/components", 20))
if type(PCA_COMPONENTS) is not int or PCA_COMPONENTS < 2:
    raise WorkflowError("pca/components (or pca/n_pcs) must be an integer of at least 2.")


if PCA_METHOD == "plink2":
    rule plink_pca_fit:
        input:
            bed=rules.create_structure_qc_set.output.bed,
            bim=rules.create_structure_qc_set.output.bim,
            fam=rules.create_structure_qc_set.output.fam,
            unrelated=rules.king_unrelated.output.keep,
        output:
            prune_in=ws_path("structure/pca/markers.prune.in"),
            prune_out=ws_path("structure/pca/markers.prune.out"),
            prune_log=ws_path("structure/pca/markers.log"),
            eigenvec=PCA_FIT_PREFIX + ".eigenvec",
            eigenval=PCA_FIT_PREFIX + ".eigenval",
            counts=PCA_FIT_PREFIX + ".acount",
            weights=PCA_FIT_PREFIX + ".eigenvec.allele",
            score_weights=ws_path("structure/pca/projection_weights.tsv"),
            log=PCA_FIT_PREFIX + ".log",
        container:
            "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
        threads: 12
        resources:
            runtime=360,
            mem_mb=65536,
        params:
            bfile=ws_path(STRUCTURE_QC_PREFIX),
            prefix=PCA_FIT_PREFIX,
            prune=ws_path("structure/pca/markers"),
            window=lambda wc: cfg("ld_pruning/window_kb", 500),
            r2=lambda wc: cfg("ld_pruning/r2", 0.2),
            maf=lambda wc: cfg("pca/plink_min_maf", cfg("thresholds/maf", 0.01)),
            n_pcs=PCA_COMPONENTS,
            seed=lambda wc: cfg("pca/seed", 1),
        shell:
            r"""
            mkdir -p "$(dirname {output.eigenvec:q})"
            workdir=$(mktemp -d "$(dirname {output.eigenvec:q})/fit.XXXXXXXX")
            trap 'rm -rf "$workdir"' EXIT

            # Re-evaluate founder status after selecting the KING-unrelated samples.
            plink2 --bfile {params.bfile:q} --keep {input.unrelated:q} \
              --make-founders --make-bed --out "$workdir/training" \
              --threads {threads} --memory {resources.mem_mb}

            n_samples=$(wc -l < "$workdir/training.fam")
            if (( n_samples < 10 )); then
                echo "Study PCA requires at least 10 KING-unrelated samples." >&2
                exit 1
            fi
            small_set=()
            if (( n_samples < 50 )); then
                small_set=(--bad-ld --bad-freqs)
            fi
            plink2 --bfile "$workdir/training" --maf {params.maf} --mac 1 --nonfounders \
              --indep-pairwise {params.window}kb 1 {params.r2} "${{small_set[@]}}" \
              --out {params.prune:q} --threads {threads} --memory {resources.mem_mb}

            n_variants=$(wc -l < {output.prune_in:q})
            k={params.n_pcs}
            if (( k >= n_samples )); then k=$((n_samples - 1)); fi
            if (( k > n_variants )); then k=$n_variants; fi
            if (( k < 2 )); then
                echo "Fewer than two markers remain after PCA pruning." >&2
                exit 1
            fi
            mode=meanimpute
            # PLINK's randomized PCA needs at least 2*k*(k+1) markers.
            if (( n_samples > 5000 && n_variants >= 2 * k * (k + 1) )); then
                mode=approx
            fi

            plink2 --bfile "$workdir/training" --extract {output.prune_in:q} \
              --nonfounders --freq counts "${{small_set[@]}}" \
              --pca allele-wts "$k" "$mode" scols=fid vcols=chrom,ref,alt \
              --out {params.prefix:q} --seed {params.seed} \
              --threads {threads} --memory {resources.mem_mb}

            # BED is biallelic. Use one ALT weight per marker for --score.
            awk 'NR == 1 || $4 == $5' {output.weights:q} > {output.score_weights:q}
            """


    rule plink_pca_project:
        input:
            bed=rules.create_structure_qc_set.output.bed,
            bim=rules.create_structure_qc_set.output.bim,
            fam=rules.create_structure_qc_set.output.fam,
            keep=rules.structure_sample_set.output.keep,
            variants=rules.plink_pca_fit.output.prune_in,
            weights=rules.plink_pca_fit.output.score_weights,
            counts=rules.plink_pca_fit.output.counts,
        output:
            scores=PCA_SCORE_PREFIX + ".sscore",
            variants=PCA_SCORE_PREFIX + ".sscore.vars",
            log=PCA_SCORE_PREFIX + ".log",
        container:
            "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
        threads: 12
        resources:
            runtime=360,
            mem_mb=32768,
        params:
            bfile=ws_path(STRUCTURE_QC_PREFIX),
            prefix=PCA_SCORE_PREFIX,
        shell:
            r"""
            last_col=$(awk 'NR == 1 {{ print NF; exit }}' {input.weights:q})
            plink2 --bfile {params.bfile:q} --keep {input.keep:q} \
              --extract {input.variants:q} --read-freq {input.counts:q} \
              --score {input.weights:q} 2 5 header-read variance-standardize list-variants cols=fid,scoresums \
              --score-col-nums "6-$last_col" --out {params.prefix:q} \
              --threads {threads} --memory {resources.mem_mb}
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
        plink=lambda wc: (
            [rules.plink_pca_fit.output.eigenvec, rules.plink_pca_fit.output.eigenval,
             rules.plink_pca_fit.output.score_weights, rules.plink_pca_project.output.scores,
             rules.plink_pca_project.output.variants]
            if PCA_METHOD == "plink2" else []
        ),
    output:
        pcs=ws_path("review/pca/pca_initial_all_samples.tsv"),
        scores_plot=ws_path("review/pca/pca_scores.pdf"),
        loadings_plot=ws_path("review/pca/pca_loadings.pdf"),
        loadings=ws_path("review/pca/pca_loadings.tsv"),
        lrldr=(
            ws_path("review/pca/autoSVD_long_range_LD_regions.tsv")
            if PCA_METHOD == "autosvd" else []
        ),
        robust=ws_path("review/pca/robust_pca_distances.tsv"),
        suggested=ws_path("review/pca/suggested_pca_outliers.tsv"),
        retained=ws_path(
            "review/pca/autoSVD_retained_variants.txt" if PCA_METHOD == "autosvd"
            else "review/pca/plink_pca_retained_variants.txt"
        ),
        summary=ws_path("review/pca/pca_initial_summary.tsv"),
    conda:
        "../envs/r_environment.yaml"
    threads: 12
    resources:
        runtime=360,
        mem_mb=65536,
    params:
        prefix=ws_path(STRUCTURE_QC_PREFIX),
        method=PCA_METHOD,
        report_input_flag=("" if PCA_METHOD == "plink2" else "--lrldr-output"),
        n_pcs=PCA_COMPONENTS,
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
          --method {params.method:q} --bed-prefix {params.prefix:q} --unrelated-file {input.unrelated:q} \
          --keep-file {input.keep:q} --n-pcs {params.n_pcs} \
          --robust-distance-n-pcs {params.robust_pcs} \
          --robust-mad-multiplier {params.robust_mad} \
          --autosvd-max-iter {params.max_iter} \
          --autosvd-min-maf {params.min_maf} --seed {params.seed} \
          --pcs-output {output.pcs:q} --scores-plot {output.scores_plot:q} \
          --loadings-plot {output.loadings_plot:q} --loadings-output {output.loadings:q} \
          --robust-output {output.robust:q} \
          --suggested-outliers {output.suggested:q} \
          --retained-variants {output.retained:q} \
          --summary-output {output.summary:q} --threads {threads} \
          {params.report_input_flag} {input.plink:q} {output.lrldr:q}
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
        echo "PCA review required: inspect scores, loadings, robust distances, suggested outliers and any method-specific outputs. This marker records completed diagnostics, not accepted exclusions. No suggested PCA outliers are removed automatically." >&2
        touch {output.ready:q}
        """