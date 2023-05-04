# LongOmics Project

## Major Reserach Project - Master Bioinfomratics and Systems Biology, Vrije Universiteit Amsterdam


The interactive application has been developed within the Long-omics project. It provides functionalities to get insights about Nanopore Long Read samples.


Files needed for the application:

</ul>
   <li> depth_LRS_ID.tsv: TSV file with coverage information. Three columns: chromosome, genomic location and count.
   <li> merged_LRS_IS.csv: Merged VCF files (merged_LRS_ID.csv), that are the output of the Variant Calling tools in VCF format, merged with the pyhton script merge_vcf.py 
   <li> sample_info.csv provides information about all the samples included in the study.
   The row names are the samples id (LRS_id). More in detail, it requires the following columns: Chr, Start, End, Length, Gene, Method (Adaptive Sampling / Cas9 Enrichment).
  </ul>
  
  
  
Snakemake scripts:
</ul>
<li> preprocessing.smk: alignment and preprocessing of fastq files with minimap2 and bedtools.  
<li> VC.smk: variant calling with SVIM, Sniffles2, CuteSV and Pepper.
<li> vep.smk and annotsv.smk: VEP and AnnotSV to assess relevance of the variants identified by the variant callers.
  </ul>


 
  
