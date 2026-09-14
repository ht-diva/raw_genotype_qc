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
    threads: 8
    resources:
        runtime=120,
        mem_mb=16384,
    params:
        plink2=PLINK2,
        bfile=ws_path(MISSING_PREFIX),
        prefix=ws_path(AUTO_QC_PREFIX),
        maf=lambda wc: cfg("variant_qc/maf", 0.01),
        geno=lambda wc: cfg("variant_qc/geno", 0.1),
        hwe=lambda wc: hwe_args(),
    shell:
        r'''
        {params.plink2} --bfile {params.bfile} --autosome \
          --maf {params.maf} --geno {params.geno} {params.hwe} \
          --make-bed --out {params.prefix} \
          --threads {threads} --memory 16000
        '''
