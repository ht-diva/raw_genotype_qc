rule create_autosomal_qc_set:
    input:
        bed=rules.apply_sample_missingness.output.bed,
        bim=rules.apply_sample_missingness.output.bim,
        fam=rules.apply_sample_missingness.output.fam,
    output:
        bed=ws_path(AUTO_QC_PREFIX + ".bed"),
        bim=ws_path(AUTO_QC_PREFIX + ".bim"),
        fam=ws_path(AUTO_QC_PREFIX + ".fam"),
        log=ws_path(AUTO_QC_PREFIX + ".log"),
    container:
        "docker://gitlab.fht.org:5050/hds-center/containers/plink2:0e8e82d8"
    threads:
        8
    resources:
        runtime=120,
        mem_mb=16384,
    params:
        bfile = ws_path(MISSING_PREFIX),
        prefix = ws_path(AUTO_QC_PREFIX),
        maf = lambda wc: cfg("variant_qc/maf"),
        hwe = lambda wc: cfg("variant_qc/hwe"),
    shell:
        r"""
        plink2 \
          --bfile {params.bfile} \
          --autosome \
          --maf {params.maf} \
          --geno {params.geno} \
          {params.hwe} \
          --make-bed \
          --out {params.prefix} \
          --threads {threads} \
          --memory {resources.mem_mb}
        """