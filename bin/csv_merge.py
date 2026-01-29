#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import glob
import os
import argparse
import numpy as np
import sys
from pathlib import Path



parser = argparse.ArgumentParser(description='Merge_vcf')

parser.add_argument('-v1', '--vcf1', dest='vcf_path1', type=str,
                        help='Sample name',nargs='?')
parser.add_argument('-v2', '--vcf2', dest='vcf_path2', type=str,
                        help='Sample name',nargs='?', const=1, default=1)
parser.add_argument('-v3', '--vcf3', dest='vcf_path3', type=str,
                        help='Sample name',nargs='?', const=1, default=1)
parser.add_argument('-v4', '--vcf4', dest='vcf_path4', type=str,
                        help='Sample name',nargs='?', const=1, default=1)
parser.add_argument('-o', '--output', dest='output_dir', type=str, default='.',
                    help='Répertoire de sortie')


args = parser.parse_args()
files =args.vcf_path1, args.vcf_path2, args.vcf_path3, args.vcf_path4
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

empty_dir = output_dir.parents[1] / "merge_intron_empty"
empty_dir.mkdir(parents=True, exist_ok=True)

samp_files = list(files)
sample_name = []
for i in samp_files:
    i = "_".join(i.split("/")[-1].split("_")[0:-1])
    if i not in sample_name:
        sample_name.append(i)

Sample_name = sample_name[0]

# joining files with concat and read_csv
Merge_df = pd.concat(map(pd.read_csv, files), ignore_index=True)
Merge_df = Merge_df.reset_index(drop=True)

df = Merge_df.replace(r'^\s*$', np.nan, regex=True)

Merge_df['Source'].str.strip()
####################

Merge_df['VarCal'] = Merge_df.groupby(["#CHROM","POS",'AA_change'])['Source'].transform(
                                              lambda x: ','.join(x))
# drop duplicate data
Merge_df = Merge_df.drop_duplicates()




################
# add confidence and Avg VAF column


CDS = Merge_df[Merge_df['Annotation'].isin(["missense_variant","synonymous_variant"])]
if CDS.empty:
    merge_empty = empty_dir / f"{Sample_name}_merge.csv"
    with open(merge_empty, "w") as fb:
        fb.write(
            "HGVS.c,POS,ALT,CHROM,AA_change,Confidence,AVG_VAF,AVG_COV,"
            "REF,Sample_name,VARTYPE,Annotation,HGVS.p,VarCal,codon\n"
        )
        fb.close()
        #sys.exit()

else:
    Merge_df_2 = CDS.groupby(["#CHROM","POS",'AA_change']).agg({'AA_change': 'count', 'VAF':'mean','DP' : 'mean'})
    Merge_df_2.columns = ['Confidence', 'AVG_VAF', 'AVG_COV']
    Merge_df_2 = Merge_df_2.reset_index()
#drop the

    """
    Merge_df_3 = CDS.drop(['Source','QUAL','Unnamed: 0','AD','DP','AF','DP4','VAF','Annotation_Impact'], axis=1)
    """
    cols_to_drop = ['Source', 'QUAL', 'Unnamed: 0', 'AD', 'DP', 'AF', 'DP4', 'VAF', 'Annotation_Impact']
    # on ne garde que les colonnes vraiment présentes
    cols_present = [c for c in cols_to_drop if c in CDS.columns]
    Merge_df_3 = CDS.drop(cols_present, axis=1)
    
    Merge_df_3 = Merge_df_3.drop_duplicates()
    
    result_1 = pd.merge(Merge_df_2, Merge_df_3,how='left',  on=["#CHROM","POS",'AA_change'])
    result_1 = result_1.rename({'#CHROM' : 'CHROM'}, axis=1)

################################
# drop duplicates while caling mnps and snps
#add new confidence

    result_1.AVG_VAF = result_1.AVG_VAF.round(2)

    df1 = result_1[result_1.duplicated(subset=['Sample_name','CHROM','POS','AVG_VAF'],keep=False)]

    df_snps = df1.loc[df1["VARTYPE"] == 'SNP']

    res = pd.merge(result_1,df_snps, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)
    res = res.replace(r'^\s*$', np.nan, regex=True)
    res['codon'] = res.AA_change.str.extract('(\d+)').astype(str)

    df2 = res[res.duplicated(subset=['Sample_name','CHROM','codon'], keep=False)]


## add new confiodence and varcal
    df3 = df2.groupby(['Sample_name','codon'])['VarCal'].agg(','.join).reset_index()

    df3['VarCal'] = df3['VarCal'].astype(str).str.split(',').apply(set).str.join(',')

    df3['Confidence_2'] =  df3["VarCal"].apply(lambda x: len(x.split(',')))
    df_merge = pd.merge(df2, df3, on=['Sample_name','codon'])

    df_merge = df_merge.drop(columns=['VarCal_x','Confidence'])
    df_merge = df_merge.rename(columns = {'VarCal_y': 'VarCal', 'Confidence_2':'Confidence' })


    
    #df2_snps = df2.loc[df2["VARTYPE"].str.contains('SNP') ]
    df2_snps = df2.loc[df2["VARTYPE"] =='SNP' ]

    res_1 = pd.merge(res, df2_snps, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)



    res_1.set_index(['HGVS.c','POS', 'ALT'],inplace=True)

    df_merge = df_merge.drop_duplicates(subset=['HGVS.c','POS','ALT'])


    res_1.update(df_merge.set_index(['HGVS.c','POS','ALT']))
    #res_1.reset_index()
    res_1.reset_index(inplace=True)
    

    #res_1.to_csv(Sample_name+"_merge"+'.csv')
    # merge output
    merge_output = os.path.join(args.output_dir, Sample_name + "_merge.csv")
    res_1.to_csv(merge_output, index=False)





###########Introns file###########


Introns = Merge_df[~Merge_df['Annotation'].isin(["missense_variant","synonymous_variant",'stop_gained'])]
#Introns = Merge_df.loc[Merge_df['AA_change'].isnull | or Merge_df['HGVS.p'].isnull()]


if Introns.empty:
    introns_empty = empty_dir / f"{Sample_name}_introns.csv"
    with open(introns_empty, "w") as fb:
        fb.write(
            "#CHROM,POS,HGVS.c,Confidence,AVG_VAF,AVG_COV,"
            "REF,ALT,Sample_name,VARTYPE,Annotation,VarCal\n"
        )
        fb.close()
        #sys.exit()
else:
        Introns['VarCal'] = Introns.groupby(["#CHROM","POS",'HGVS.c'])['Source'].transform(lambda x: ', '.join(x))

        Introns =Introns.drop_duplicates()
        Introns_1 = Introns.groupby(["#CHROM","POS",'HGVS.c']).agg({'HGVS.c': 'count', 'VAF':'mean','DP' : 'mean'})
        Introns_1.columns = ['Confidence', 'AVG_VAF', 'AVG_COV']
        Introns_1 = Introns_1.reset_index()
        Introns = Introns.drop(['Source','QUAL','Unnamed: 0','AD','DP','AF','DP4','VAF','Annotation_Impact','AA_change','HGVS.p'], axis=1)
        Introns = Introns.drop_duplicates()

        result_2 = pd.merge(Introns_1, Introns ,how='left',  on=["#CHROM","POS",'HGVS.c'])
        result_2 = result_2[result_2.Confidence > 2]

        #result_2.to_csv(Sample_name+"_introns"+'.csv')
        # introns output
        introns_output = os.path.join(args.output_dir, Sample_name + "_introns.csv") 
        result_2.to_csv(introns_output, index=False)



        print(f"[OK] Fichier généré : {merge_output}, {introns_output}")


"""
#!/usr/bin/env python3


import pandas as pd
import os
import argparse
import numpy as np
from pathlib import Path

parser = argparse.ArgumentParser(description='Fusionne les fichiers CSV de différents callers')
parser.add_argument('-v1', '--vcf1', dest='vcf_path1', type=str, required=True, help='Chemin du CSV Samtools')
parser.add_argument('-v2', '--vcf2', dest='vcf_path2', type=str, required=True, help='Chemin du CSV FreeBayes')
parser.add_argument('-v3', '--vcf3', dest='vcf_path3', type=str, required=True, help='Chemin du CSV HaplotypeCaller')
parser.add_argument('-v4', '--vcf4', dest='vcf_path4', type=str, required=True, help='Chemin du CSV VarDict')
parser.add_argument('-o', '--output', dest='output_dir', type=str, default='.', help='Répertoire de sortie')

args = parser.parse_args()

# --- Préparation des chemins ---
files = [args.vcf_path1, args.vcf_path2, args.vcf_path3, args.vcf_path4]
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

# --- Nom de l’échantillon ---
sample_name = "_".join(Path(files[0]).stem.split("_")[0:-1])
if not sample_name:
    sample_name = "sample"

merge_file = output_dir / f"{sample_name}_merge.csv"
introns_file = output_dir / f"{sample_name}_introns.csv"

try:
    # --- Lecture et concaténation ---
    Merge_df = pd.concat(map(pd.read_csv, files), ignore_index=True)
    Merge_df = Merge_df.replace(r'^\s*$', np.nan, regex=True)
    Merge_df['Source'] = Merge_df['Source'].astype(str).str.strip()

    # --- VarCal ---
    Merge_df['VarCal'] = Merge_df.groupby(["#CHROM", "POS", "AA_change"])['Source'].transform(lambda x: ','.join(x))
    Merge_df = Merge_df.drop_duplicates()

    # --- CDS : Missense et Synonymous ---
    CDS = Merge_df[Merge_df['Annotation'].isin(["missense_variant", "synonymous_variant"])]

    if CDS.empty:
        # Créer fichier vide avec en-tête
        with open(merge_file, "w") as fb:
            fb.write("HGVS.c,POS,ALT,CHROM,AA_change,Confidence,AVG_VAF,AVG_COV,REF,Sample_name,VARTYPE,Annotation,HGVS.p,VarCal,codon\n")
    else:
        Merge_df_2 = CDS.groupby(["#CHROM", "POS", 'AA_change']).agg({'AA_change': 'count', 'VAF': 'mean', 'DP': 'mean'})
        Merge_df_2.columns = ['Confidence', 'AVG_VAF', 'AVG_COV']
        Merge_df_2 = Merge_df_2.reset_index()

        Merge_df_3 = CDS.drop(['Source', 'QUAL', 'Unnamed: 0', 'AD', 'DP', 'AF', 'DP4', 'VAF', 'Annotation_Impact'], axis=1, errors='ignore')
        Merge_df_3 = Merge_df_3.drop_duplicates()

        result_1 = pd.merge(Merge_df_2, Merge_df_3, how='left', on=["#CHROM", "POS", 'AA_change'])
        result_1 = result_1.rename({'#CHROM': 'CHROM'}, axis=1)
        result_1['AVG_VAF'] = result_1['AVG_VAF'].round(2)
        result_1['codon'] = result_1['AA_change'].str.extract('(\d+)')


        # Sauvegarde du merge principal
        result_1.to_csv(merge_file, index=False)

    # --- Introns ---
    Introns = Merge_df[~Merge_df['Annotation'].isin(["missense_variant", "synonymous_variant", 'stop_gained'])]

    if Introns.empty:
        with open(introns_file, "w") as fb:
            fb.write("#CHROM,POS,HGVS.c,Confidence,AVG_VAF,AVG_COV,REF,ALT,Sample_name,VARTYPE,Annotation,VarCal\n")
    else:
        Introns['VarCal'] = Introns.groupby(["#CHROM", "POS", 'HGVS.c'])['Source'].transform(lambda x: ','.join(x))
        Introns = Introns.drop_duplicates()

        Introns_1 = Introns.groupby(["#CHROM", "POS", 'HGVS.c']).agg({'HGVS.c': 'count', 'VAF': 'mean', 'DP': 'mean'})
        Introns_1.columns = ['Confidence', 'AVG_VAF', 'AVG_COV']
        Introns_1 = Introns_1.reset_index()

        Introns = Introns.drop(
            ['Source', 'QUAL', 'Unnamed: 0', 'AD', 'DP', 'AF', 'DP4', 'VAF', 'Annotation_Impact', 'AA_change', 'HGVS.p'],
            axis=1, errors='ignore'
        ).drop_duplicates()

        result_2 = pd.merge(Introns_1, Introns, how='left', on=["#CHROM", "POS", 'HGVS.c'])
        result_2 = result_2[result_2['Confidence'] > 2]

        result_2.to_csv(introns_file, index=False)

    print(f"[OK] {sample_name} traité avec succès → {merge_file.name}, {introns_file.name}")

except Exception as e:
    print(f"[ERREUR] {sample_name} : {e}")
    raise SystemExit(1)
"""