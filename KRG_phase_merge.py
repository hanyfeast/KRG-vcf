#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import scipy as sp
import re
import matplotlib.pyplot as plt
import seaborn as sns


cwd = "KorRefGenome/"
phase1 = "variants0622"
phase2 = "variants1100"


def align_indel_left(seq1,seq2, pos):
    if len(seq1) > len(seq2):
        a = seq1
        b = seq2
    else:
        a = seq2
        b = seq1

    for i in range(1,len(b)+1):
        if a[i*-1] != b[i*-1]:
            i = i - 1
            for j in range(1,len(b)-i+1):
                if (a[j] != b[j]) | (j>len(b)-i):
                    break
            break
    
    if len(seq2) > len(seq1):
        return pd.Series([pos+j-1, b[j-1], a[j-1:i*-1]], index=['new_pos','ref.Allele','alt.Allele'])
    else:
        return pd.Series([pos+j-1, a[j-1:i*-1], b[j-1]], index=['new_pos','ref.Allele','alt.Allele'])

for i in range(1,25):
    if i == 23:
        i = "X"
    elif i == 24:
        i = "Y"

    ##### phase1 ######
    phase1_chr_cmm = pd.read_csv(cwd + phase1 + "_cmm_chr" + str(i) + ".txt",sep="\t",header=0)
    phase1_chr_rare = pd.read_csv(cwd + phase1 + "_rare_chr" + str(i) + ".txt",sep="\t",header=0)
    phase1_chr_indel = pd.read_csv(cwd + phase1 + "_indel_chr" + str(i) + ".txt",sep="\t",header=0)
    
    phase2_chr_cmm = pd.read_csv(cwd + phase2 + "_cmm_chr" + str(i) + ".txt",sep="\t",header=0)
    phase2_chr_rare = pd.read_csv(cwd + phase2 + "_rare_chr" + str(i) + ".txt",sep="\t",header=0)
    phase2_chr_indel = pd.read_csv(cwd + phase2 + "_indel_chr" + str(i) + ".txt",sep="\t",header=0)
    
    #phase1 indel parsing
    phase1_chr_indel['alt2'] = phase1_chr_indel['{alt}'].str.rstrip(to_strip=";")
    phase1_chr_indel['alt.cnt2'] = phase1_chr_indel['{alt.cnt}'].str.rstrip(to_strip=";")
    phase1_chr_indel['alt.freq2'] = phase1_chr_indel['{alt.freq}'].str.rstrip(to_strip=";")

    a1 = phase1_chr_indel['alt2'].str.split(";",expand=True).stack().reset_index(level=1,drop=True)
    a2 = phase1_chr_indel['alt.cnt2'].str.split(";",expand=True).stack().reset_index(level=1,drop=True)
    a3 = phase1_chr_indel['alt.freq2'].str.split(";",expand=True).stack().reset_index(level=1,drop=True)
    phase1_chr_indel_split = pd.concat([a1,a2,a3],axis=1, keys=['alt','alt.cnt','alt.freq'])
    
    phase1_chr_indel2 = phase1_chr_indel.merge(phase1_chr_indel_split, left_index=True, right_index=True, how='left')
    
    phase1_chr_indel2['ref.Allele'] = ""
    phase1_chr_indel2['alt.Allele'] = ""
    phase1_chr_indel2['new_pos'] = phase1_chr_indel2['pos']
    
    phase1_chr_indel2_insertion = phase1_chr_indel2[phase1_chr_indel2['reference'].str.len()<phase1_chr_indel2['alt'].str.len()]
    phase1_chr_indel2_deletion = phase1_chr_indel2[phase1_chr_indel2['reference'].str.len()>phase1_chr_indel2['alt'].str.len()]
    phase1_chr_indel2_delins = phase1_chr_indel2[phase1_chr_indel2['reference'].str.len()==phase1_chr_indel2['alt'].str.len()]

    print("p1_ins",sum(phase1_chr_indel2_insertion['ref.Allele']==""))
    sub1 = ((phase1_chr_indel2_insertion.apply(lambda x: x['reference'][1:] == x['alt'][(len(x['reference'])-1)*-1:],axis=1) | 
        (phase1_chr_indel2_insertion['reference'].str.len()==1)) & (phase1_chr_indel2_insertion['alt'].str[0] == phase1_chr_indel2_insertion['reference'].str[0]))
    phase1_chr_indel2_insertion.loc[sub1, "alt.Allele"] = phase1_chr_indel2_insertion[sub1].apply(lambda x: x['alt'][:len(x['alt']) - len(x['reference']) + 1],axis=1)
    phase1_chr_indel2_insertion.loc[sub1, 'ref.Allele'] = phase1_chr_indel2_insertion[sub1]['reference'].str[0]
    print("p1_ins",sum(phase1_chr_indel2_insertion['ref.Allele']==""))
    
    print("p1_del",sum(phase1_chr_indel2_deletion['ref.Allele']==""))
    sub1 = ((phase1_chr_indel2_deletion.apply(lambda x: x['alt'][1:] == x['reference'][(len(x['alt'])-1)*-1:],axis=1) | 
        (phase1_chr_indel2_deletion['alt'].str.len()==1)) & (phase1_chr_indel2_deletion['reference'].str[0] == phase1_chr_indel2_deletion['alt'].str[0]))
    phase1_chr_indel2_deletion.loc[sub1, "ref.Allele"] = phase1_chr_indel2_deletion[sub1].apply(lambda x: x['reference'][:len(x['reference']) - len(x['alt']) + 1],axis=1)
    phase1_chr_indel2_deletion.loc[sub1, "alt.Allele"] = phase1_chr_indel2_deletion[sub1]['alt'].str[0]
    print("p1_del",sum(phase1_chr_indel2_deletion['ref.Allele']==""))
    
    print("p1_delins",sum(phase1_chr_indel2_delins['ref.Allele']==""))
    
    sub1 = (phase1_chr_indel2_delins['reference'].str[1:] == phase1_chr_indel2_delins['alt'].str[1:]) | (phase1_chr_indel2_delins['reference'].str.len()==1)
    phase1_chr_indel2_delins.loc[sub1,'ref.Allele'] = phase1_chr_indel2_delins[sub1]['reference'].str[0]
    phase1_chr_indel2_delins.loc[sub1,'alt.Allele'] = phase1_chr_indel2_delins[sub1]['alt'].str[0]
    
    phase1_chr_indel2_delins['code'] = phase1_chr_indel2_delins['chr'] + "_" + phase1_chr_indel2_delins['pos'].astype(str) + "_" + phase1_chr_indel2_delins['ref.Allele'] + "_" + phase1_chr_indel2_delins['alt.Allele']
    phase1_chr_indel2_delins = phase1_chr_indel2_delins[
        ~(phase1_chr_indel2_delins['code'].isin(phase1_chr_cmm['chr'] + "_" + phase1_chr_cmm['pos'].astype(str) + "_" + phase1_chr_cmm['ref.Allele'] + "_" + phase1_chr_cmm['alt.Allele']) | 
          phase1_chr_indel2_delins['code'].isin(phase1_chr_rare['chr'] + "_" + phase1_chr_rare['pos'].astype(str) + "_" + phase1_chr_rare['ref.Allele'] + "_" + phase1_chr_rare['alt.Allele']))
    ]
    print("p1_delins",sum(phase1_chr_indel2_delins['ref.Allele']==""))

    phase1_chr_indel3 = pd.concat([phase1_chr_indel2_insertion, phase1_chr_indel2_deletion, phase1_chr_indel2_delins])
    phase1_chr_indel3 = phase1_chr_indel3.sort_index()
    
    phase1_chr_cmm['ref.freq'] = 1 - phase1_chr_cmm['alt.freq']
    phase1_chr_rare['ref.freq'] = 1 - phase1_chr_rare['alt.freq']
    
    phase1_chr_cmm['new_pos'] = phase1_chr_cmm['pos']
    phase1_chr_rare['new_pos'] = phase1_chr_rare['pos']
    
    phase1_all = pd.concat([
        phase1_chr_cmm[['chr','new_pos','ref.Allele','alt.Allele','ref.freq','alt.freq']],
        phase1_chr_rare[['chr','new_pos','ref.Allele','alt.Allele','ref.freq','alt.freq']],
        phase1_chr_indel3[['chr','new_pos','ref.Allele','alt.Allele','ref.freq','alt.freq']]
    ], ignore_index=True)
    phase1_all.columns = ['chr','pos','ref.Allele','alt.Allele','ref.freq','alt.freq']
    
    phase1_all_single = phase1_all[['chr','pos','ref.Allele','alt.Allele','ref.freq','alt.freq']].drop_duplicates(subset=['chr','pos','ref.Allele','alt.Allele'],keep=False)
    
    phase1_all_multi = phase1_all[phase1_all[['chr','pos','ref.Allele','alt.Allele','ref.freq','alt.freq']].duplicated(subset=['chr','pos','ref.Allele','alt.Allele'],keep=False)]
    phase1_all_multi = phase1_all_multi[['chr','pos','ref.Allele','alt.Allele','ref.freq','alt.freq']].sort_values(by=['pos','alt.freq'],ascending=[True,False])
    phase1_all_multi_wid = phase1_all_multi[~phase1_all_multi.duplicated(subset=['chr','pos','ref.Allele','alt.Allele'],keep='first')].merge(
        phase1_all_multi[~phase1_all_multi.duplicated(subset=['chr','pos','ref.Allele','alt.Allele'],keep='last')],how="left",on=['chr','pos','ref.Allele','alt.Allele'])
    phase1_all_multi_wid['ref.freq'] = phase1_all_multi_wid['ref.freq_x'] - phase1_all_multi_wid['alt.freq_y'].astype(float)
    phase1_all_multi_wid['alt.freq'] = phase1_all_multi_wid['alt.freq_x'].astype(float) + phase1_all_multi_wid['alt.freq_y'].astype(float)
    phase1_all_final = pd.concat([phase1_all_single,phase1_all_multi_wid[['chr','pos','ref.Allele','alt.Allele','ref.freq','alt.freq']]])
    
    phase1_all_final['ref.cnt'] = (phase1_all_final['ref.freq'].astype(float) * 622 * 2).round(0)
    phase1_all_final['alt.cnt'] = (phase1_all_final['alt.freq'].astype(float) * 622 * 2).round(0)
    
    
    #########################################################
    ######################## phase2 #########################
    #########################################################
    
    
    phase2_chr_cmm['alt_info'] = phase2_chr_cmm['{alt:freq}'].str.rstrip(to_strip=", ")
    phase2_chr_cmm_split = phase2_chr_cmm['alt_info'].str.split(", ",expand=True).stack().reset_index(level=1,drop=True).to_frame('alt2')
    phase2_chr_cmm2 = phase2_chr_cmm.merge(phase2_chr_cmm_split, left_index=True, right_index=True, how='left')
    phase2_chr_cmm2[['alt.Allele','alt.freq']] = phase2_chr_cmm2['alt2'].str.split(":",expand=True)

    phase2_chr_rare['alt_info'] = phase2_chr_rare['{alt:freq}'].str.rstrip(to_strip=", ")
    phase2_chr_rare_split = phase2_chr_rare['alt_info'].str.split(", ",expand=True).stack().reset_index(level=1,drop=True).to_frame('alt2')
    phase2_chr_rare2 = phase2_chr_rare.merge(phase2_chr_rare_split, left_index=True, right_index=True, how='left')
    phase2_chr_rare2[['alt.Allele','alt.freq']] = phase2_chr_rare2['alt2'].str.split(":",expand=True)

    phase2_chr_cmm2['new_pos'] = phase2_chr_cmm2['pos']
    phase2_chr_rare2['new_pos'] = phase2_chr_rare2['pos']
    phase2_chr_cmm2['ref.Allele'] = phase2_chr_cmm2['reference']
    phase2_chr_rare2['ref.Allele'] = phase2_chr_rare2['reference']
    
    phase2_chr_indel['alt_info'] = phase2_chr_indel['{alt:freq}'].str.rstrip(to_strip=", ")
    phase2_chr_indel_split = phase2_chr_indel['alt_info'].str.split(", ",expand=True).stack().reset_index(level=1,drop=True).to_frame('alt2')
    phase2_chr_indel2 = phase2_chr_indel.merge(phase2_chr_indel_split, left_index=True, right_index=True, how='left')
    phase2_chr_indel2[['alt','alt.freq']] = phase2_chr_indel2['alt2'].str.split(":",expand=True)
    
    phase2_chr_indel2['new_pos'] = phase2_chr_indel2['pos']
    phase2_chr_indel2['ref.Allele'] = ""
    phase2_chr_indel2['alt.Allele'] = ""
    
    phase2_chr_indel2_insertion = phase2_chr_indel2[phase2_chr_indel2['reference'].str.len()<phase2_chr_indel2['alt'].str.len()]
    phase2_chr_indel2_deletion = phase2_chr_indel2[phase2_chr_indel2['reference'].str.len()>phase2_chr_indel2['alt'].str.len()]
    phase2_chr_indel2_delins = phase2_chr_indel2[phase2_chr_indel2['reference'].str.len()==phase2_chr_indel2['alt'].str.len()]
    
    print("p2_ins",sum(phase2_chr_indel2_insertion['ref.Allele']==""))
    #ref length 1일때 또는 ref length >1 일때 alt와 왼쪽 1, 오른쪽 나머지가 양끝으로 일치할때
    sub1 = ((phase2_chr_indel2_insertion.apply(lambda x: x['reference'][1:] == x['alt'][(len(x['reference'])-1)*-1:],axis=1) | 
        (phase2_chr_indel2_insertion['reference'].str.len()==1)) & (phase2_chr_indel2_insertion['alt'].str[0] == phase2_chr_indel2_insertion['reference'].str[0]))

    phase2_chr_indel2_insertion.loc[sub1, "ref.Allele"] = phase2_chr_indel2_insertion[sub1]['reference'].str[0]
    phase2_chr_indel2_insertion.loc[sub1, "alt.Allele"] = phase2_chr_indel2_insertion[sub1].apply(lambda x: x['alt'][:len(x['alt']) - len(x['reference']) + 1],axis=1)
    print("p2_ins",sum(phase2_chr_indel2_insertion['ref.Allele']==""))
    
    #reference가 alt외 왼쪽부터 일치할때
    sub1 = phase2_chr_indel2_insertion.apply(lambda x: x['alt'][:len(x['reference'])] == x['reference'],axis=1) & (phase2_chr_indel2_insertion['ref.Allele'] == "")
    phase2_chr_indel2_insertion.loc[sub1, "ref.Allele"] = phase2_chr_indel2_insertion[sub1]['reference'].str[-1]
    phase2_chr_indel2_insertion.loc[sub1, "alt.Allele"] = phase2_chr_indel2_insertion[sub1].apply(lambda x: x['alt'][len(x['reference'])-1:],axis=1)
    phase2_chr_indel2_insertion.loc[sub1, "new_pos"] = phase2_chr_indel2_insertion[sub1]['pos'] + phase2_chr_indel2_insertion[sub1]['reference'].str.len() - phase2_chr_indel2_insertion[sub1]['ref.Allele'].str.len()
    print("p2_ins",sum(phase2_chr_indel2_insertion['ref.Allele']==""))
    
    #reference가 alt외 좌우로 일치할때
    sub1 = (phase2_chr_indel2_insertion['ref.Allele'] == "")
    phase2_chr_indel2_insertion.loc[sub1,['new_pos','ref.Allele','alt.Allele']] = phase2_chr_indel2_insertion[sub1].apply(lambda x:align_indel_left(x['reference'],x['alt'],x['pos']),axis=1)
    print("p2_ins",sum(phase2_chr_indel2_insertion['ref.Allele']==""))

    print("p2_del",sum(phase2_chr_indel2_deletion['ref.Allele']==""))
    sub1 = ((phase2_chr_indel2_deletion.apply(lambda x: x['alt'][1:] == x['reference'][(len(x['alt'])-1)*-1:],axis=1) | 
        (phase2_chr_indel2_deletion['alt'].str.len()==1)) & (phase2_chr_indel2_deletion['reference'].str[0] == phase2_chr_indel2_deletion['alt'].str[0]))
    phase2_chr_indel2_deletion.loc[sub1, "ref.Allele"] = phase2_chr_indel2_deletion[sub1].apply(lambda x: x['reference'][:len(x['reference']) - len(x['alt']) + 1],axis=1)
    phase2_chr_indel2_deletion.loc[sub1, "alt.Allele"] = phase2_chr_indel2_deletion[sub1]['alt'].str[0]
    print("p2_del",sum(phase2_chr_indel2_deletion['ref.Allele']==""))
    sub1 = phase2_chr_indel2_deletion.apply(lambda x: x['reference'][:len(x['alt'])] == x['alt'],axis=1) & (phase2_chr_indel2_deletion['ref.Allele'] == "")
    phase2_chr_indel2_deletion.loc[sub1, "ref.Allele"] = phase2_chr_indel2_deletion[sub1].apply(lambda x: x['reference'][len(x['alt'])-1:],axis=1)
    phase2_chr_indel2_deletion.loc[sub1, "alt.Allele"] = phase2_chr_indel2_deletion[sub1]['alt'].str[-1]
    phase2_chr_indel2_deletion.loc[sub1, "new_pos"] = phase2_chr_indel2_deletion[sub1]['pos'] + phase2_chr_indel2_deletion[sub1]['reference'].str.len() - phase2_chr_indel2_deletion[sub1]['ref.Allele'].str.len()
    print("p2_del",sum(phase2_chr_indel2_deletion['ref.Allele']==""))
    sub1 = (phase2_chr_indel2_deletion['ref.Allele'] == "")
    phase2_chr_indel2_deletion.loc[sub1,['new_pos','ref.Allele','alt.Allele']] = phase2_chr_indel2_deletion[sub1].apply(lambda x:align_indel_left(x['reference'],x['alt'],x['pos']),axis=1)
    print("p2_del",sum(phase2_chr_indel2_deletion['ref.Allele']==""))
    
    #delins는 SNV형식이거나 왼쪽 1개만 다르고 다 똑같음
    print("p2_delins",sum(phase2_chr_indel2_delins['ref.Allele']==""))
    sub1 = (phase2_chr_indel2_delins['reference'].str[1:] == phase2_chr_indel2_delins['alt'].str[1:]) | (phase2_chr_indel2_delins['reference'].str.len()==1)
    phase2_chr_indel2_delins.loc[sub1,'ref.Allele'] = phase2_chr_indel2_delins[sub1]['reference'].str[0]
    phase2_chr_indel2_delins.loc[sub1,'alt.Allele'] = phase2_chr_indel2_delins[sub1]['alt'].str[0]
    
    #이미 있는 SNP의 경우 제거
    phase2_chr_indel2_delins['code'] = phase2_chr_indel2_delins['chr'] + "_" + phase2_chr_indel2_delins['pos'].astype(str) + "_" + phase2_chr_indel2_delins['ref.Allele'] + "_" + phase2_chr_indel2_delins['alt.Allele']
    phase2_chr_indel2_delins = phase2_chr_indel2_delins[
        ~(phase2_chr_indel2_delins['code'].isin(phase2_chr_cmm2['chr'] + "_" + phase2_chr_cmm2['pos'].astype(str) + "_" + phase2_chr_cmm2['ref.Allele'] + "_" + phase2_chr_cmm2['alt.Allele']) | 
          phase2_chr_indel2_delins['code'].isin(phase2_chr_rare2['chr'] + "_" + phase2_chr_rare2['pos'].astype(str) + "_" + phase2_chr_rare2['ref.Allele'] + "_" + phase2_chr_rare2['alt.Allele']))
    ] 
    print("p2_delins",sum(phase2_chr_indel2_delins['ref.Allele']==""))
    
    phase2_chr_indel3 = pd.concat([phase2_chr_indel2_insertion, phase2_chr_indel2_deletion, phase2_chr_indel2_delins])
    phase2_chr_indel3 = phase2_chr_indel3.sort_index()
    
    phase2_all = pd.concat([
        phase2_chr_cmm2[['chr','new_pos','ref.Allele','alt.Allele','ref.freq','alt.freq']],
        phase2_chr_rare2[['chr','new_pos','ref.Allele','alt.Allele','ref.freq','alt.freq']],
        phase2_chr_indel3[['chr','new_pos','ref.Allele','alt.Allele','ref.freq','alt.freq']]
    ], ignore_index=True)
    phase2_all.columns = ['chr','pos','ref.Allele','alt.Allele','ref.freq','alt.freq']
    
    phase2_all_single = phase2_all[['chr','pos','ref.Allele','alt.Allele','ref.freq','alt.freq']].drop_duplicates(subset=['chr','pos','ref.Allele','alt.Allele'],keep=False)
    
    phase2_all_multi = phase2_all[phase2_all[['chr','pos','ref.Allele','alt.Allele','ref.freq','alt.freq']].duplicated(subset=['chr','pos','ref.Allele','alt.Allele'],keep=False)]
    phase2_all_multi = phase2_all_multi[['chr','pos','ref.Allele','alt.Allele','ref.freq','alt.freq']].sort_values(by=['pos','alt.freq'],ascending=[True,False])
    phase2_all_multi_wid = phase2_all_multi[~phase2_all_multi.duplicated(subset=['chr','pos','ref.Allele','alt.Allele'],keep='first')].merge(
        phase2_all_multi[~phase2_all_multi.duplicated(subset=['chr','pos','ref.Allele','alt.Allele'],keep='last')],how="left",on=['chr','pos','ref.Allele','alt.Allele'])
    phase2_all_multi_wid['ref.freq'] = phase2_all_multi_wid['ref.freq_x'] - phase2_all_multi_wid['alt.freq_y'].astype(float)
    phase2_all_multi_wid['alt.freq'] = phase2_all_multi_wid['alt.freq_x'].astype(float) + phase2_all_multi_wid['alt.freq_y'].astype(float)
    phase2_all_final = pd.concat([phase2_all_single,phase2_all_multi_wid[['chr','pos','ref.Allele','alt.Allele','ref.freq','alt.freq']]])
    
    phase2_all_final['ref.cnt'] = (phase2_all_final['ref.freq'].astype(float) * 1100 * 2).round(0)
    phase2_all_final['alt.cnt'] = (phase2_all_final['alt.freq'].astype(float) * 1100 * 2).round(0)
    
    phase1_2_all = pd.merge(phase1_all_final, phase2_all_final, how='outer', on=['chr','pos','ref.Allele','alt.Allele'], suffixes=('_p1','_p2'))
    phase1_2_all = phase1_2_all.fillna(0)
    
    print(i, len(phase1_all_final), len(phase2_all_final), len(phase1_2_all))
    
    phase1_2_all['ref.cnt'] = phase1_2_all['ref.cnt_p1'] + phase1_2_all['ref.cnt_p2']
    phase1_2_all['alt.cnt'] = phase1_2_all['alt.cnt_p1'] + phase1_2_all['alt.cnt_p2']
    phase1_2_all['ref.freq'] = phase1_2_all['ref.cnt']/(622+1100)/2
    phase1_2_all['alt.freq'] = phase1_2_all['alt.cnt']/(622+1100)/2
    phase1_2_all.loc[(phase1_2_all['ref.cnt_p2'] == 0) & (phase1_2_all['alt.cnt_p2'] == 0),'ref.freq'] = phase1_2_all[(phase1_2_all['ref.cnt_p2'] == 0) & (phase1_2_all['alt.cnt_p2'] == 0)]['ref.cnt_p1']/622/2
    phase1_2_all.loc[(phase1_2_all['ref.cnt_p2'] == 0) & (phase1_2_all['alt.cnt_p2'] == 0),'alt.freq'] = phase1_2_all[(phase1_2_all['ref.cnt_p2'] == 0) & (phase1_2_all['alt.cnt_p2'] == 0)]['alt.cnt_p1']/622/2
    phase1_2_all.loc[(phase1_2_all['ref.cnt_p1'] == 0) & (phase1_2_all['alt.cnt_p1'] == 0),'ref.freq'] = phase1_2_all[(phase1_2_all['ref.cnt_p1'] == 0) & (phase1_2_all['alt.cnt_p1'] == 0)]['ref.cnt_p2']/1100/2
    phase1_2_all.loc[(phase1_2_all['ref.cnt_p1'] == 0) & (phase1_2_all['alt.cnt_p1'] == 0),'alt.freq'] = phase1_2_all[(phase1_2_all['ref.cnt_p1'] == 0) & (phase1_2_all['alt.cnt_p1'] == 0)]['alt.cnt_p2']/1100/2
    
    phase1_2_all = phase1_2_all.sort_values(by=['pos','alt.freq'],ascending=[True,False])
    phase1_2_all[['pos','ref.cnt_p1','alt.cnt_p1','ref.cnt_p2','alt.cnt_p2','ref.cnt','alt.cnt']] = phase1_2_all[['pos','ref.cnt_p1','alt.cnt_p1','ref.cnt_p2','alt.cnt_p2','ref.cnt','alt.cnt']].astype(int)
    phase1_2_all.to_csv("KRG_merge_song_chr" + str(i) + "_sorted.csv",sep="\t",index=False,header=True)

vcfHeader1 = """##fileformat=VCFv4.1
##fileDate=20191026
##source=KRG_merge_table_to_vcf
##contig=<ID="""
vcfHeader2 = """>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=RN1,Number=1,Type=Integer,Description="Reference count in Phase1">
##INFO=<ID=AN1,Number=1,Type=Integer,Description="Alt count in Phase1">
##INFO=<ID=RN2,Number=1,Type=Integer,Description="Reference count in Phase2">
##INFO=<ID=AN2,Number=1,Type=Integer,Description="Alt count in Phase2">
##INFO=<ID=RN0,Number=1,Type=Integer,Description="Reference count in P1+P2">
##INFO=<ID=AN0,Number=1,Type=Integer,Description="Alt count in P1+P2">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""

for i in range(1,25):
    if i == 23:
        i = "X"
    elif i == 24:
        i = "Y"
        
    variant_table = pd.read_csv("KRG_merge_song_chr" + str(i) + "_sorted.csv",sep="\t",header=0)
    vcfFile = open("KRG_merge_song_chr" + str(i) + "_sorted.vcf", "w")
    vcfFile.write(vcfHeader1 + str(i) + vcfHeader2)
    for index, row in variant_table.iterrows():
        vcfLine = "\t".join(map(str,[
            row['chr'][3:],row['pos'],".",row['ref.Allele'],row['alt.Allele'],".","PASS",
            "NS=" + str(int((row['ref.cnt'] + row['alt.cnt'])/2)) + 
            ";AF=" + str(row['alt.freq']) + 
            ";RN1=" + str(int(row['ref.cnt_p1'])) + 
            ";AN1=" + str(int(row['alt.cnt_p1'])) + 
            ";RN2=" + str(int(row['ref.cnt_p2'])) + 
            ";AN2=" + str(int(row['alt.cnt_p2'])) + 
            ";RN0=" + str(int(row['ref.cnt'])) + 
            ";AN0=" + str(int(row['alt.cnt']))
                  ]))
        vcfFile.write(vcfLine + "\n")
    vcfFile.close()
    print(i, len(variant_table))
