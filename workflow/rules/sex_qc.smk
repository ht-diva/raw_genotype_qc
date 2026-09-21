rule split_par_for_sex_qc:
    input:
        bed=rules.apply_heterozygosity_exclusions.output.bed,
        bim=rules.apply_heterozygosity_exclusions.output.bim,
        fam=rules.apply_heterozygosity_exclusions.output.fam,
    output:
        bed=ws_path("sex_qc/split_par/genotype.bed"),
        bim=ws_path("sex_qc/split_par/genotype.bim"),
        fam=ws_path("sex_qc/split_par/genotype.fam"),
        log=ws_path("sex_qc/split_par/genotype.log"),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=90,
        mem_mb=12000,
    params:
        bfile=ws_path("heterozygosity/genotype"),
        prefix=ws_path("sex_qc/split_par/genotype"),
        build=lambda wc: cfg("genome_build", "b38"),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.bed}")"

        plink2 \
            --bfile "{params.bfile}" \
            --split-par "{params.build}" \
            --make-bed \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory 11000
        """

rule create_x_qc_markers:
    input:
        bed=rules.split_par_for_sex_qc.output.bed,
        bim=rules.split_par_for_sex_qc.output.bim,
        fam=rules.split_par_for_sex_qc.output.fam,
    output:
        bed=ws_path("sex_qc/x_qc/genotype.bed"),
        bim=ws_path("sex_qc/x_qc/genotype.bim"),
        fam=ws_path("sex_qc/x_qc/genotype.fam"),
        log=ws_path("sex_qc/x_qc/genotype.log"),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=90,
        mem_mb=12000,
    params:
        bfile=ws_path("sex_qc/split_par/genotype"),
        prefix=ws_path("sex_qc/x_qc/genotype"),
        maf=lambda wc: cfg("thresholds/maf"),
        geno=lambda wc: cfg("thresholds/variant_missingness"),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.bed}")"

        plink2 \
            --bfile "{params.bfile}" \
            --chr X \
            --maf "{params.maf}" \
            --geno "{params.geno}" \
            --make-bed \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory 11000
        """

rule check_sex_candidate_thresholds:
    input:
        bed=rules.create_x_qc_markers.output.bed,
        bim=rules.create_x_qc_markers.output.bim,
        fam=rules.create_x_qc_markers.output.fam,
    output:
        sexcheck=ws_path(
            "sex_qc/sex_candidate.sexcheck"
        ),
        log=ws_path(
            "sex_qc/sex_candidate.log"
        ),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=60,
        mem_mb=8000,
    params:
        bfile=ws_path(
            "sex_qc/x_qc/genotype"
        ),
        prefix=ws_path(
            "sex_qc/sex_candidate"
        ),
        female_max=lambda wc: cfg(
            "thresholds/female_max_f"
        ),
        male_min=lambda wc: cfg(
            "thresholds/male_min_f"
        ),
    shell:
        r"""
        set -euo pipefail

        plink2 \
            --bfile "{params.bfile}" \
            --check-sex \
                max-female-xf="{params.female_max}" \
                min-male-xf="{params.male_min}" \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory 7500
        """

rule sex_qc_review:
    input:
        sexcheck=(
            rules.check_sex_candidate_thresholds.output.sexcheck
        ),
        fam=rules.create_x_qc_markers.output.fam,
        unmatched_genotypes=rules.prepare_reported_sex.output.unmatched_genotypes,
    output:
        plot=ws_path(
            "review/sex/sex_F_distribution.pdf"
        ),
        table=ws_path(
            "review/sex/sex_candidate_classification.tsv"
        ),
        candidate_thresholds=ws_path(
            "review/sex/candidate_sex_thresholds.yaml"
        ),
        summary=ws_path(
            "review/sex/sex_qc.summary.tsv"
        ),
    conda:
        "../envs/r_environment.yaml"
    params:
        female_max=lambda wc: cfg(
            "thresholds/female_max_f"
        ),
        male_min=lambda wc: cfg(
            "thresholds/male_min_f"
        ),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.plot}")"

        Rscript workflow/scripts/sex_qc_review.R \
            --sexcheck "{input.sexcheck}" \
            --fam "{input.fam}" \
            --unmatched-genotypes "{input.unmatched_genotypes}" \
            --female-max-f "{params.female_max}" \
            --male-min-f "{params.male_min}" \
            --plot "{output.plot}" \
            --table "{output.table}" \
            --candidate-thresholds \
                "{output.candidate_thresholds}" \
            --summary "{output.summary}"
        """

rule apply_approved_sex_thresholds:
    input:
        sexcheck=(
            rules.check_sex_candidate_thresholds.output.sexcheck
        ),
        fam=rules.create_x_qc_markers.output.fam,
        unmatched_genotypes=rules.prepare_reported_sex.output.unmatched_genotypes,
        threshold=lambda wc: cfg(
            "accepted_sex_thresholds",
            "config/accepted_sex_thresholds.yaml",
        ),
    output:
        exclusions=ws_path(
            "sex_qc/automatic_sex_exclusions.tsv"
        ),
        classification=ws_path(
            "sex_qc/sex_classification.tsv"
        ),
        summary=ws_path(
            "sex_qc/accepted_sex_qc.summary.tsv"
        ),
    conda:
        "../envs/r_environment.yaml"
    shell:
        r"""
        set -euo pipefail

        Rscript workflow/scripts/apply_sex_thresholds.R \
            --sexcheck "{input.sexcheck}" \
            --fam "{input.fam}" \
            --unmatched-genotypes "{input.unmatched_genotypes}" \
            --thresholds "{input.threshold}" \
            --exclusions "{output.exclusions}" \
            --classification "{output.classification}" \
            --summary "{output.summary}"
        """
