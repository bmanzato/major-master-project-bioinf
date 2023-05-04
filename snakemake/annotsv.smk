
list_LRS = []



OUTPUT_PATH = "/home/benedettamanzato/longomics_app/variant_effect_analysis/annotsv"



rule all:
    input: [OUTPUT_PATH + "/annotsv_{LRS_id}.tsv".format(LRS_id = lrs) for lrs in list_LRS]


rule depth:
    input: "/home/benedettamanzato/longomics_app/vc/{LRS_id}.vcf"
    output: "/home/benedettamanzato/longomics_app/variant_effect_analysis/annotsv/annotsv_{LRS_id}.tsv"
    shell: "./AnnotSV -SVinputFile {input} -outputFile {output}"

        
        