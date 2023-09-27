import subprocess
import pandas
from collections import OrderedDict


TARGET_REPORT = [
    str(QC_FP/'reports'/'preprocess_summary.tsv'),
    str(QC_FP/'reports'/'fastqc_quality.tsv')
]


TARGET_FASTQC = expand(
    str(QC_FP/'reports'/'{sample}_{rp}_fastqc'/'fastqc_data.txt'),
    sample=Samples.keys(), rp=Pairs)

rule all_reports:
    input:
        TARGET_REPORT


rule fastqc_report:
    """ make fastqc reports """
    input:
        files = expand(
            str(QC_FP/'04_fastqc'/'{sample}_{rp}_fastqc/fastqc_data.txt'),
            sample=Samples.keys(),rp=Pairs)
    output:
        str(QC_FP/'reports'/'fastqc_quality.tsv')
    group:
        "report"
    script:
        "../scripts/fastqc_report.py"


rule per_host_reads:
    input:
        expand(str(QC_FP/'03_decontam'/'ids'/'{{host}}'/'{sample}'),
            sample=Samples.keys()),
    output:
        str(QC_FP/'log'/'decontam'/'{host}_counts.tsv')
    params:
        decontam_fp = str(QC_FP/'03_decontam'/'ids'/'{host}')
    group:
        "report"
    shell:
        """
        set +o pipefail
        wc -l {params.decontam_fp}/* | grep -v "total" | awk '{{print $2, $1}}' | \
            awk -v OFS="\t" '{{sub(".*/", "", $1); print $0}}' > {output}
        sed -i '1i sample\t{wildcards.host}' {output}
        """


rule total_host_reads:
    input:
        expand(str(QC_FP/'03_decontam'/'hostreads'/'{sample}'), sample=Samples.keys())
    output:
        str(QC_FP/'log'/'decontam'/'host_counts.tsv')
    params:
        hostids_fp = str(QC_FP/'03_decontam'/'hostreads')
    group:
        "report"
    shell:
        """
        set +o pipefail
        wc -l {params.hostids_fp}/* | grep -v "total" | awk '{{print $2, $1}}' | \
            awk -v OFS="\t" '{{sub(".*/", "", $1); print $0}}' > {output}
        sed -i '1i sample\thost' {output}
        """


rule preprocess_report:
    """Combines the information from multiple preprocessing steps"""
    input:
        trim_files = expand(
            str(QC_FP/'log'/'trimmomatic'/'{sample}.out'),
            sample=sorted(Samples.keys())),
        bbduk_files = expand(
            str(QC_FP/'log'/'bbduk'/'{sample}.log'),
            sample=sorted(Samples.keys())),
        decontam_files = expand(
            str(QC_FP/'log'/'decontam'/'{column}_counts.tsv'),
            column=list(HostGenomes.keys())+['host'])
    output:
        str(QC_FP/'reports'/'preprocess_summary.tsv')
    group:
        "report"
    resources:
        walltime_hr=1
    script:
        "../scripts/preprocess_report.py"
