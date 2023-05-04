
import pandas as pd
import allel
import numpy as np


list_LRS = ["LRS_42", "LRS_43", "LRS_44", "LRS_46", "LRS_47", "LRS_49", "LRS_50","LRS_51",
            "LRS_53", "LRS_54","LRS_55","LRS_56", "LRS_62", "LRS_70", "LRS_71"]


for sample in list_LRS:

    # cutesv
    cutesv = allel.vcf_to_dataframe('vc/'+sample+'_cutesv.vcf',fields='*')
    if not str(cutesv) == 'None':
        cutesv = cutesv.loc[cutesv['PRECISE'] == True]
        cutesv['Tool'] = 'cuteSV'
        cutesv['QUAL'] = 'Precise'

    # sniffles2
    sniffles2 = allel.vcf_to_dataframe('vc/'+sample+'_sniffles.vcf',fields='*')
    sniffles2 = sniffles2.loc[sniffles2['QUAL'] >= 20]
    if not str(sniffles2) == 'None':
        sniffles2['Tool'] = 'Sniffles2'

    # svim
    svim = allel.vcf_to_dataframe('vc/'+sample+'_svim.vcf',fields='*')
    svim = svim.loc[svim['QUAL'] >= 20]
    if not str(svim) == 'None':
        svim['Tool'] = 'SVIM'


    # pepper
    pepper = allel.vcf_to_dataframe('vc/'+sample+'_pepper.vcf', fields='*')
    pepper = pepper.loc[pepper['QUAL'] > 30]
    pepper['ID'] = 'Pepper_' + sample

    for i in pepper.index:
        if pepper.loc[i,'is_snp'] == True:
            pepper.loc[i,'SVTYPE'] = 'SNP'
            pepper.loc[i,'SVLEN'] = 1
            pepper.loc[i,'END'] = pepper.loc[i,'POS']
        if pepper.loc[i, 'is_snp'] == False:
            if len(pepper.loc[i,'REF']) > len(pepper.loc[i,'ALT_1']):
                pepper.loc[i,'SVTYPE'] = 'DEL'
            else:
                pepper.loc[i,'SVTYPE'] = 'INS'
            pepper.loc[i,'SVLEN'] = abs(len(pepper.loc[i,'REF']) - len(pepper.loc[i,'ALT_1']))
            pepper.loc[i, 'END'] = int(pepper.loc[i, 'POS']+pepper.loc[i,'SVLEN'])
    pepper['QUAL'] = np.trunc(pepper['QUAL'])



    VCF = pd.concat([sniffles2, svim, pepper, cutesv], axis=0)

    VCF.to_csv('vc/merged_'+sample+'.csv')














