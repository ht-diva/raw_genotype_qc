rule compute_heterozygosity:
    input:
        bed=rules.create_autosomal_qc_set.output.bed,
        bim=rules.create_autosomal_qc_set.output.bim,
        fam=rules.create_autosomal_qc_set.output.fam,
    output:
        het=ws_path(
            "heterozygosity/heterozygosity.het"
        ),
        log=ws_path(
            "heterozygosity/heterozygosity.log"
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
            "autosomal_qc/genotype"
        ),
        prefix=ws_path(
            "heterozygosity/heterozygosity"
        ),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.het}")"

        plink2 \
            --bfile "{params.bfile}" \
            --het \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory 7500
        """


rule heterozygosity_qc:
    input:
        het=rules.compute_heterozygosity.output.het,
    output:
        plot=ws_path(
            "heterozygosity/heterozygosity.png"
        ),
        exclusions=ws_path(
            "heterozygosity/heterozygosity_exclusions.tsv"
        ),
        table=ws_path(
            "heterozygosity/heterozygosity_all_samples.tsv"
        ),
        summary=ws_path(
            "heterozygosity/heterozygosity.summary.tsv"
        ),
    conda:
        "../envs/r_environment.yaml"
    params:
        n_sd=lambda wc: cfg(
            "thresholds/heterozygosity_sd"
        ),
    shell:
        r"""
        set -euo pipefail

        Rscript workflow/scripts/heterozygosity_qc.R \
            --het "{input.het}" \
            --n-sd "{params.n_sd}" \
            --plot "{output.plot}" \
            --suggested-exclusions "{output.exclusions}" \
            --table "{output.table}" \
            --summary "{output.summary}"
        """


rule apply_heterozygosity_exclusions:
    input:
        bed=rules.apply_sample_missingness.output.bed,
        bim=rules.apply_sample_missingness.output.bim,
        fam=rules.apply_sample_missingness.output.fam,
        exclusions=rules.heterozygosity_qc.output.exclusions,
    output:
        bed=ws_path(
            "heterozygosity/genotype.bed"
        ),
        bim=ws_path(
            "heterozygosity/genotype.bim"
        ),
        fam=ws_path(
            "heterozygosity/genotype.fam"
        ),
        log=ws_path(
            "heterozygosity/genotype.log"
        ),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=90,
        mem_mb=12000,
    params:
        bfile=ws_path(
            "sample_missingness/genotype"
        ),
        prefix=ws_path(
            "heterozygosity/genotype"
        ),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.bed}")"

        plink2 \
            --bfile "{params.bfile}" \
            --remove "{input.exclusions}" \
            --make-bed \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory 11000
        """