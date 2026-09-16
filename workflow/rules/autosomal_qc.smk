rule create_autosomal_qc_set:
    input:
        bed=rules.apply_sample_missingness.output.bed,
        bim=rules.apply_sample_missingness.output.bim,
        fam=rules.apply_sample_missingness.output.fam,
    output:
        bed=ws_path(
            "autosomal_qc/genotype.bed"
        ),
        bim=ws_path(
            "autosomal_qc/genotype.bim"
        ),
        fam=ws_path(
            "autosomal_qc/genotype.fam"
        ),
        log=ws_path(
            "autosomal_qc/genotype.log"
        ),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=120,
        mem_mb=16384,
    params:
        bfile=ws_path(
            "sample_missingness/genotype"
        ),
        prefix=ws_path(
            "autosomal_qc/genotype"
        ),
        maf=lambda wc: cfg(
            "thresholds/maf"
        ),
        geno=lambda wc: cfg(
            "thresholds/variant_missingness"
        ),
        hwe=lambda wc: cfg(
            "thresholds/hwe"
        ),
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.bed}")"

        plink2 \
            --bfile "{params.bfile}" \
            --autosome \
            --maf "{params.maf}" \
            --geno "{params.geno}" \
            --hwe "{params.hwe}" \
            --make-bed \
            --out "{params.prefix}" \
            --threads {threads} \
            --memory {resources.mem_mb}
        """