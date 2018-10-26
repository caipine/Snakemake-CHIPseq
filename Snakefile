shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

localrules: all
# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs

# load cluster config file
CLUSTER = json.load(open(config['CLUSTER_JSON']))
FILES = json.load(open(config['SAMPLES_JSON']))

import csv
import os

SAMPLES = sorted(FILES.keys())

## list all samples by sample_name and sample_type
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        MARK_SAMPLES.append(sample + "_" + sample_type)


# which sample_type is used as control for calling peaks: e.g. Input, IgG...
CONTROL = config["control"]
CONTROLS = [sample for sample in MARK_SAMPLES if CONTROL in sample]
CASES = [sample for sample in MARK_SAMPLES if CONTROL not in sample]

## multiple samples may use the same control input/IgG files
CONTROLS_UNIQUE = list(set(CONTROLS))


## list BAM files
CONTROL_BAM = expand("03aln/{sample}.sorted.bam", sample=CONTROLS_UNIQUE)
CASE_BAM = expand("03aln/{sample}.sorted.bam", sample=CASES)

## peaks and bigwigs
ALL_PEAKS = []
ALL_inputSubtract_BIGWIG = []
ALL_SUPER = []

for case in CASES:
    sample = "_".join(case.split("_")[0:-1])
    control = sample + "_" + CONTROL
    if control in CONTROLS:
        ALL_PEAKS.append("08peak_macs1/{}_vs_{}_macs1_peaks.bed".format(case, control))
        ALL_PEAKS.append("08peak_macs1/{}_vs_{}_macs1_nomodel_peaks.bed".format(case, control))
        ALL_PEAKS.append("09peak_macs2/{}_vs_{}_macs2_peaks.xls".format(case, control))
        ALL_inputSubtract_BIGWIG.append("06bigwig_inputSubtract/{}_subtract_{}.bw".format(case, control))
        ALL_SUPER.append("11superEnhancer/{}_vs_{}-super/".format(case, control))

ALL_SAMPLES = CASES + CONTROLS_UNIQUE
ALL_BAM     = CONTROL_BAM + CASE_BAM
ALL_DOWNSAMPLE_BAM = expand("04aln_downsample/{sample}-downsample.sorted.bam", sample = ALL_SAMPLES)
ALL_FASTQ   = expand("01seq/{sample}.fastq", sample = ALL_SAMPLES)
ALL_FASTQC  = expand("02fqc/{sample}_fastqc.zip", sample = ALL_SAMPLES)
ALL_INDEX = expand("03aln/{sample}.sorted.bam.bai", sample = ALL_SAMPLES)
ALL_DOWNSAMPLE_INDEX = expand("04aln_downsample/{sample}-downsample.sorted.bam.bai", sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand("03aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES)
ALL_PHATOM = expand("05phantompeakqual/{sample}_phantom.txt", sample = ALL_SAMPLES)
ALL_BIGWIG = expand("07bigwig/{sample}.bw", sample = ALL_SAMPLES)
ALL_QC = ["10multiQC/multiQC_log.html"]



TARGETS = []
TARGETS.extend(ALL_FASTQC)

## sometimes if if you have TF ChIP-seq data, do not include it to chromHMM, or you want
## only a subset of the histone marks be included in the chromHMM call

if config["chromHMM"]:
    HISTONE_INCLUDED = config["histone_for_chromHMM"].split(" ")
    HISTONE_CASES = [sample for sample in MARK_SAMPLES if sample.split("_")[-1] in HISTONE_INCLUDED ]
    ALL_BED = expand("12bed/{sample}.bed", sample = HISTONE_CASES + CONTROLS)
    CHROMHMM = ["13chromHMM/MYOUTPUT", "13chromHMM/binarizedData"]
    CHROMHMM_TABLE = ["12bed/cellmarkfiletable.txt"]

localrules: all
rule all:
     input: TARGETS


## get a list of fastq.gz files for the same mark, same sample
def get_fastq(wildcards):
    sample = "_".join(wildcards.sample.split("_")[0:-1])
    mark = wildcards.sample.split("_")[-1]
    return FILES[sample][mark]

## now only for single-end ChIPseq,
rule merge_fastqs:
    input: get_fastq
    output: "01seq/{sample}.fastq"
    log: "00log/{sample}_unzip"
    threads: CLUSTER["merge_fastqs"]["cpu"]
    params: jobname = "{sample}"
    message: "merging fastqs gunzip -c {input} > {output}"
    shell: "gunzip -c {input} > {output} 2> {log}"


rule fastqc:
    input: "01seq/{sample}.fastq"
    output: "02fqc/{sample}_fastqc.zip", "02fqc/{sample}_fastqc.html"
    log:    "00log/{sample}_fastqc"
    threads: CLUSTER["fastqc"]["cpu"]
    params : jobname = "{sample}"
    message: "fastqc {input}: {threads}"
    shell:
        """
        fastqc -o 02fqc -f fastq --noextract {input} 2> {log}
        """

######################################################
# get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and dupliated reads by samblaster -r
# samblaster should run before samtools sort
rule align:
    input: "01seq/{sample}.fastq"
    output: "03aln/{sample}.sorted.bam", "00log/{sample}.align"
    threads: CLUSTER["align"]["cpu"]
    params:
            bowtie = "--chunkmbs 320 -m 1 --best -p 16 ",
            jobname = "{sample}"
    message: "aligning {input}: {threads} threads"
    log:
        bowtie = "00log/{sample}.align",
        markdup = "00log/{sample}.markdup"
    shell:
        """
        bowtie {params.bowtie} {config[idx_bt1]} -q {input} -S 2> {log.bowtie} \
        | samblaster --removeDups \
	| samtools view -Sb -F 4 - \
	| samtools sort -m 2G -@ 5 -T {output[0]}.tmp -o {output[0]} 2> {log.markdup}

        """



