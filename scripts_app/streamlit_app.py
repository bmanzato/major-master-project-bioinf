
import os
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import streamlit as st
from PIL import Image


st.set_page_config(
    page_title='Long-Omics',
    layout='centered',
    initial_sidebar_state='expanded',
    menu_items={
            'About': "Long-Omics\n\n" +
            "Questions? 🙋 " +
            "Feedback? 📝 " +
            "Found a bug? \n\n" +
            "Please contact us at contact@lizard.bio"
    })


image = Image.open('bioLizard_logo.png')
st.image(image, width=140)


st.sidebar.write('### Long-Omics')
explanation = """
<p>
    This interactive application has been developed within the Long-omics project.
</p>
<p>
    The Coverage Plot shows the coverage of the target region of the samples belonging to the chosen phenotype.
    Structural variants and genes spanning the target regions can be plotted as well.
</p>
<p>
</ul>
<hr>
"""

st.sidebar.markdown(explanation, unsafe_allow_html=True)

st.title('Long-Omics')


def coverage_plot(target_start, target_end):

    #st.text('Target region of this plot is: ' + str(target_start) + ':' + str(target_end))

    # slice the count data with the preferred target region
    data_new = data.loc[target_start : target_end,:]
    target_lenght = target_end - target_start

    # make a subset of the SV data so that SV present in the target region will be shown in the table
    VCF = data_vcf[data_vcf['CHROM'].isin([c,chr])]
    VCF = VCF.loc[VCF['POS'] >= target_start-500]
    VCF = VCF.loc[VCF['POS'] <= target_end+500]
    VCF = VCF.loc[VCF['SVTYPE'] != 'BND']

    VCF_display = VCF
    VCF_display.index = VCF_display['ID']
    VCF_display = VCF_display[['CHROM','POS','REF','ALT_1','QUAL','SVLEN','SVTYPE','END']]

    Overview_Variants = st.checkbox('Overview of the variants identified for sample ' + sample, False) 
    if Overview_Variants:
            st.dataframe(VCF_display, use_container_width=False)

    # Coverage Plot #############################
    if target_lenght >= 600000: n = 150
    elif 450000 < target_lenght < 600000: n = 100
    elif 6000 < target_lenght < 450000: n = 16
    elif 950 < target_lenght < 6000: n = 5
    else: n = 1

     # average the counts, n is dependent on the lenght of the target region
    data_chr_cov = data_new.groupby(np.arange(len(data_new)) // n).mean()

    fig = px.bar(data_chr_cov, x='Location', y='Count', #title='Coverage Plot Chromosome ' + str(c) + ' - ' + sample,
                 color='Count', labels={'Location', 'Count'}, color_continuous_scale = "darkmint")

    lower = target_start - target_lenght / 20
    upper = target_end + target_lenght / 20

    MAX = max(data['Count'])
    ylim = MAX + MAX / 8

    fig.update_xaxes(range=[lower, upper], showgrid=False)
    fig.update_yaxes(range=[0, ylim], visible=True)

    ymin = 0.2

    # Gene Plot #################################
    target_genes['colors'] = target_genes['Gene name'].apply(lambda x: 'DarkRed' if x == sample_info.loc[sample, 'Gene'] else 'Grey')
    w2 = st.checkbox('Show ' + sample + ' genes above the coverage plot'   , False)
    if w2:
        for index, row in target_genes.iterrows():
            fig.add_trace(go.Scatter(x=[row['Gene start (bp)'], row['Gene end (bp)']], y=[ymin, ymin],
                                 mode='lines', line=dict(color=row.colors, width=5), name=row['Gene name'],
                                 opacity=1))


    # SV Plot #################################
    # this plots the SVs above the coverage plot
    if not VCF.empty:
        VCF['colors'] = VCF['SVTYPE'].apply(lambda x: 'Blue' if x == 'DUP' else ('Red' if x == 'INS' else ('Orange' if x == 'INV' else 'Green')))
        for index, row in VCF.iterrows():
            fig.add_trace(go.Scatter(x=[row['POS'], row['POS']+abs(row['SVLEN'])], y=[MAX + MAX/10, MAX + MAX/10],
                                     mode='lines',line=dict(color=row.colors, width=5), name=row['ID'], opacity=1))


    # start and end of target region (the red vertical lines)
    fig.add_shape(type="line", x0=target_start, y0=0,x1=target_start, y1=ylim, line=dict(color="DarkRed", width=1))
    fig.add_shape(type="line", x0=target_end, y0=0,x1=target_end, y1=ylim, line=dict(color="DarkRed", width=1))


    fig.update_layout(autosize=False, width=900, height=300, margin=dict(l=40, r=20, b=50, t=30, pad=4),
                      showlegend=False, xaxis=dict(showgrid=False),yaxis=dict(showgrid=False))


    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)', 'paper_bgcolor': 'rgba(0, 0, 0, 0)', })

    st.plotly_chart(fig)


    avg_targetreg = sum(data['Count']) / len(data['Count'])
    if sample_info_viz.loc[phenotype, 'Gene'] in list(target_genes.index):
        data_gene1 = data.loc[data['Location'] >= genes.loc[sample_info.loc[sample,'Gene']]['Gene start (bp)']]
        data_gene1 = data_gene1.loc[data['Location'] <= genes.loc[sample_info.loc[sample, 'Gene']]['Gene end (bp)']]
        avg_gene1 = sum(data_gene1['Count']) / len(data_gene1['Count'])
        st.text('The average coverage of the target gene ' + sample_info.loc[sample, 'Gene'] + ' is ' + str(
            round(avg_gene1, 2)))


    zeros = (data['Count'] == 0).sum()

    st.text('Sample ' + sample + ' has an average coverage of ' + str(round(avg_targetreg,2)) +
            '. \nNumber of positions with zero counts is ' + str(zeros) + ' (' + str(round(zeros/len(data['Count']),2)*100) + '%)')

    N = VCF.shape[0]
    ins = VCF.loc[VCF['SVTYPE'] == 'INS'].shape[0]
    dell = VCF.loc[VCF['SVTYPE'] == 'DEL'].shape[0]
    inv = VCF.loc[VCF['SVTYPE'] == 'INV'].shape[0]

    st.text(str(N) + ' SVs are found for sample ' + sample)
    if N != 0:
        st.text(str(ins) + ' (' + str(round(ins / N,2) * 100) + '%) of them are Insertions. \n' +
                str(dell) + ' (' + str(round(dell / N, 2) * 100) + '%) of them are Deletions. \n' +
                str(inv) + ' (' + str(round(inv / N, 2) * 100) + '%) of them are Inversions. \n')


    # Zoom in plot for each SV:
    # Detailed plots ############################
    if not VCF.empty:
        choice = list(VCF.index)
        svs = st.multiselect(
            'Would you like to look further into any SV for sample ' + sample + ' ?', choice)

        vcf_choice = VCF.loc[svs, :]
        st.text('SV chosen for further analysis:')
        st.dataframe(vcf_choice[['CHROM','POS','REF','ALT_1','QUAL','SVLEN','SVTYPE','END']], use_container_width=False)

        # add colors to vcf_choice:
        vcf_choice['colors'] = vcf_choice['SVTYPE'].apply(lambda x: 'Blue' if x == 'DUP' else ('Red' if x == 'INS' else ('Orange' if x == 'INV' else 'Green')))

        for index, row in vcf_choice.iterrows():

            if row['SVLEN']+2000 >= 360000:
                n = 120
            elif 360000 < row['SVLEN']+2000 < 160000:
                n = 80
            elif 50000 < row['SVLEN']+2000 < 160000:
                n = 60
            elif 2200 < row['SVLEN']+2000 < 50000:
                n = 25
            else:
                n = 15

            # average the counts, n is dependent on the lenght of the target region
            data_sv = data_new.groupby(np.arange(len(data_new)) // n).mean()


            sv_fig = px.bar(data_sv, x='Location', y='Count', title='ZOOM IN Coverage Plot Chromosome ' + str(c) +
                                                                     ' - ' + sample + ' SV: ' + row['ID'],
                        color='Count',labels={'Location', 'Count'}, color_continuous_scale = "darkmint")

            # lower and upper limits
            lower = row['POS'] - 8000
            upper = row['POS'] + abs(row['SVLEN']) + 8000

            MAX = max(data_sv['Count'])
            ylim = MAX + MAX / 8

            sv_fig.update_xaxes(range=[lower, upper], showgrid=False)
            sv_fig.update_yaxes(range=[0, ylim], visible=True)

            ymin = 0.2

            # Gene Plot #################################
            for index, row in target_genes.iterrows():
                fig.add_trace(go.Scatter(x=[row['Gene start (bp)'], row['Gene end (bp)']], y = [ymin, ymin],
                                         mode ='lines', line = dict(color=row.colors, width=5), name = row['Gene name'], opacity = 1))
                fig.add_trace(go.Scatter(x=[row['Gene start (bp)'], row['Gene end (bp)']], y = [ymin, ymin],
                                         mode ='markers + text', name ='Markers and Text',text = row['Gene name'], textposition ='bottom center'))

            # SV Plot ###################################
            for index, row in vcf_choice.iterrows():
                sv_fig.add_trace(
                    go.Scatter(x=[row['POS'], row['POS'] + abs(row['SVLEN'])], y=[MAX + MAX / 10, MAX + MAX / 10],
                           mode='lines', line=dict(color=row.colors, width=5), text='', name=row['ID'], opacity=1))

            sv_fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)', 'paper_bgcolor': 'rgba(0, 0, 0, 0)', })
            sv_fig.update_layout(autosize=False, width=900, height=300, margin=dict(l=40, r=20, b=100, t=30, pad=4),
                               showlegend=False, xaxis=dict(showgrid=False))

            st.plotly_chart(sv_fig)




def gene_plot(goi,target_start,target_end): # genes of interest; list

    VCF = data_vcf[data_vcf['CHROM'].isin([c, chr])]
    VCF = VCF.loc[VCF['POS'] >= target_start - 500]
    VCF = VCF.loc[VCF['POS'] <= target_end + 500]

    for i in goi:

        gl = genes.loc[i, 'Gene end (bp)'] - genes.loc[i, 'Gene start (bp)']
        gene_lenght = 40 * gl

        # data_gene = data[~((data['Location'] > genes.loc[i, 'Gene start (bp)'] - 20*gl) & (data['Location'] < genes.loc[i, 'Gene end (bp)'] + 20*gl))]
        data_gene = data.loc[genes.loc[i, 'Gene start (bp)'] - 20 * gl : genes.loc[i, 'Gene end (bp)'] + 20 * gl]
        gene_fig = px.bar(data_gene, x='Location', y='Count', title='ZOOM IN Gene ' + i, color='Count',
                          labels={'Location', 'Count'}, height=400, color_continuous_scale = "darkmint")

        # lower and upper limits
        lower = genes.loc[i, 'Gene start (bp)'] - 0.2*gl
        upper = genes.loc[i, 'Gene end (bp)'] + 0.2*gl

        MAX = max(data_gene['Count'])
        ylim = MAX + MAX / 8

        gene_fig.update_xaxes(range=[lower, upper], showgrid=False)
        gene_fig.update_yaxes(range=[0, ylim], visible=True)

        # Gene
        ymin = 0.2
        gene_fig.add_trace(go.Scatter(x=[genes.loc[i, 'Gene start (bp)'], genes.loc[i, 'Gene end (bp)']],
                                      y=[0.2, ymin],mode='lines', line=dict(color='DarkRed', width=5), opacity=1))

        # SV
        for index, row in VCF.iterrows():
            gene_fig.add_trace(go.Scatter(x=[row['POS'], row['POS'] + abs(row['SVLEN'])],
                                          y=[MAX + MAX / 10, MAX + MAX / 10], mode='lines',
                                          #line=dict(color=row.colors, width=5),
                                          text='', name=row['ID'], opacity=1))

        gene_fig.update_layout(autosize=False, width=900, height=300, margin=dict(l=40, r=20, b=100, t=30, pad=4),
                               showlegend=False, xaxis=dict(showgrid=False),yaxis=dict(showgrid=False))

        gene_fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)', 'paper_bgcolor': 'rgba(0, 0, 0, 0)', })

        st.plotly_chart(gene_fig)


uploaded_file = st.file_uploader('Upload the sample data file:', help=('Please provide the sample data as a .csv or '
    'Excel file (.xlsx or .xls). For each sample, it requires the following columns: '
    'Chr, Start, End, Length, Gene, Method (Adaptive Sampling / Cas9 Enrichment).'))

if uploaded_file is not None:
    file_name, file_extension = os.path.splitext(uploaded_file.name)
    if file_extension == '.csv':
        try:
            sample_info = pd.read_csv(uploaded_file, sep=',', header=0, index_col=0)
            st.write('*Sample Data upload...* ✅')
        except:
            st.error(f'Could not read the input file {uploaded_file.name}')
    elif file_extension == '.xlsx' or file_extension == '.xls':
        try:
            sample_info = pd.read_excel(uploaded_file, header=0, index_col=0)
            st.write('*Sample Data upload...* ✅')
        except Exception as e:
            st.error(f'Could not read the input file {uploaded_file.name}')

    sample_info_viz = sample_info.drop_duplicates('Phenotype')
    sample_info_viz.index = sample_info_viz['Phenotype']
    sample_info_viz.drop('Phenotype',inplace=True,axis=1)
    #st.dataframe(sample_info_viz, use_container_width=False)
    Overview_Samples = st.checkbox('Overview of samples', False)
    if Overview_Samples:
        st.dataframe(sample_info_viz, use_container_width=False)

    # Choose a phenotype
    phenotype = st.selectbox(
        'Select a phenotype to be analyzed:',
        ('BPES (non-coding deletion)', 'BPES (Poly-Ala expansion)', 'NF1', 'Huntington Disease',
         'Fragile X Syndrome', 'Steinert Myotonic Dystrophy','Polycystische nierziekte','DSD'))

    sample_info = sample_info.loc[sample_info['Phenotype'] == phenotype]

    chr = 'chr' + sample_info.iloc[0, 0]
    c = sample_info.iloc[0, 0]

    samples_list = list(sample_info.index)

    st.text('The target Region of ' + sample_info.iloc[0, 6] + ' is: Chromosome ' + str(c) + ':' +
            str(sample_info.iloc[0, 1]) + '-' + str(sample_info.iloc[0, 2]) +
            '.\nLenght is ' + str(sample_info.iloc[0, 3]) + 'b.')

    
    mart = pd.read_csv('mart_export.txt',sep='\t',index_col=0,header=0)
    mart = mart.loc[mart['Chromosome/scaffold name'] == sample_info_viz.loc[phenotype,'Chr']]
    mart = mart.dropna()
    mart = mart.drop('Transcript stable ID',axis=1)

    # Find genes spanning the target region
    gene_indexes = []
    for index, row in mart.iterrows():
        if sample_info.iloc[0, 1] < row['Gene start (bp)'] < sample_info.iloc[0, 2]:
            gene_indexes.append(index)
        if sample_info.iloc[0, 1] < row['Gene end (bp)']:
            if sample_info.iloc[0, 2] > row['Gene end (bp)']:#+100:
                if index not in gene_indexes:
                    gene_indexes.append(index)

    target_genes = mart.loc[gene_indexes, :]
    target_genes = target_genes.drop_duplicates()

    target_genes['Lenght'] = target_genes['Gene end (bp)'] - target_genes['Gene start (bp)']
    target_genes = target_genes.rename(columns={'Chromosome/scaffold name': 'Chr'})

    genes = target_genes
    genes.index = genes['Gene name'] # HGNC symbol
    genes = genes.drop('Gene name', axis=1)
    Overview_Genes = st.checkbox('Overview of the genes spanning the target region', False) 
    if Overview_Genes:
        st.dataframe(genes,use_container_width=False)

    ##
    st.text(phenotype + ' phenotype has ' + str(len(sample_info)) + ' samples.')
    
    choice_sample=list(samples_list)
    sample = st.selectbox('For which sample would you like to show the coverage plot? ', choice_sample)

    
    data = pd.read_csv('depth/depth_' + sample + '.tsv', sep='\t', header=None)
    data.rename(columns={0: 'Chr', 1: 'Location', 2: 'Count'}, inplace=True)
    data.index = data['Location']
    data = data.dropna()
    #data['avg'] = sum(data['Count']) / len(data['Count'])
    
    data_vcf = pd.read_csv('vc/merged_' + sample + '.csv', sep=',', header=0, index_col=0)
    # Coverage Plots

    coverage_plot(sample_info.loc[sample, 'Start'], sample_info.loc[sample, 'End'])


    if sample_info_viz.loc[phenotype, 'Gene'] in list(target_genes.index):
        st.subheader('Genes plot')
        Display_Plot_Gene = st.checkbox('Display Gene plot', False)
        if Display_Plot_Gene:
            gene_plot([sample_info_viz.loc[phenotype, 'Gene']],sample_info.loc[sample, 'Start'],sample_info.loc[sample, 'End'])


 