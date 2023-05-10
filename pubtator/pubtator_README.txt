PubMed_query.txt is the link used to find publications on PubMed related to BPES.


Output of that step is in bpes_syndrome_queries.txt. The queries are splitted in blocks of 100 IDs, as the PubTator API only allows a maximum of 100 PubMed IDs per query.


The output of PubTator are Title, Abstract and Table with Genes and Mutations mentioned in the paper. The results are concatenated into a single file: pubtator_results_bpes_syndrome.txt.


Quick analysis of results in bpes_variant_analysis.ipynb. The scripts organizes in a table form the PubTator results.
Additionally, in the notebook mutations_df can be saved and translated to the preferred mutation name (g. / p. etc) with Transvar (https://bioinformatics.mdanderson.org/transvar/)