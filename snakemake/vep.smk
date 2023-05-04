list_LRS = []


OUTPUT_PATH = "/home/benedettamanzato/longomics_app/variant_effect_analysis/vep"


rule all:
    input: [OUTPUT_PATH + "/vep_{LRS_id}.tsv".format(LRS_id = lrs) for lrs in LRS_tools]
    
    
rule vep:
    input: "/home/benedettamanzato/longomics_app/vc/{LRS_id}.vcf"
    output: "/home/benedettamanzato/longomics_app/variant_effect_analysis/vep/vep_{LRS_id}.tsv"
    shell: "./vep --af --appris --biotype --buffer_size 500 --check_existing --distance 5000 --mane --pubmed --regulatory --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file {input} --output_file {output} --force_overwrite --refseq --chr chr3 --format vcf"