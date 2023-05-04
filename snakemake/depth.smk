list_LRS = []


OUTPUT_PATH = "/home/benedettamanzato/longomics_app/depth"


rule all:
    input: [OUTPUT_PATH + "/depth_{LRS_id}.tsv".format(LRS_id = lrs) for lrs in list_LRS]


rule align:
    input: "/home/benedettamanzato/data_longomics/fastq/{LRS_id}.fastq.gz"
    output: "/home/benedettamanzato/data_longomics/bam_files/{LRS_id}.bam"
    shell: "minimap2 -ax map-ont ref_fastq/grch38.fasta.gz {input} | samtools view -b -F 256 -F 2048 {input} > {output}"
    

rule sort:
    input: "/home/benedettamanzato/data_longomics/bam_files/{LRS_id}.bam"
    output: "/home/benedettamanzato/data_longomics/bam_files/{LRS_id}_sorted.bam"
    shell: "samtools sort {input} > {output}"


rule index:
    input: "/home/benedettamanzato/data_longomics/bam_files/{LRS_id}_sorted.bam"
    output: "/home/benedettamanzato/data_longomics/bam_files/{LRS_id}.bam.bai"
    shell: "samtools index {input} > {output}"
    
    
rule depth:
    input: "/home/benedettamanzato/data_longomics/bam_files/{LRS_id}_sorted.bam"
    output: "/home/benedettamanzato/data_longomics/depth/depth_{LRS_id}.tsv"
    shell: "samtools depth -a -b targets.bed -Q 0 {input} > {output}"