#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 14:34:23 2021

@author: yunda_si
"""


import os
import cluster_species as cs
import pickle
import rw_a2m
import sys


def final_pair(TaxID_dict,topn):
    
    paired = []
    
    for TaxID in TaxID_dict:
        faA_list,faB_list = TaxID_dict[TaxID]

        faA_list = faA_list[:topn]
        faB_list = faB_list[:topn]
        
        for faA,faB in zip(faA_list,faB_list):             

            header = faA[-1][0].strip()+'||'+faB[-1][0].strip()
            seq = faA[-1][1].strip()+faB[-1][1].strip()
            paired.append([header,seq])  
        
    print(f'len paired: {len(paired)}')
    return paired


def main(file_dict,cov,topn):

    file = open(file_dict['tax2id'],'rb')
    tax2id = pickle.load(file)
    file.close()

    file = open(file_dict['tax_db'],'rb')
    tax_db = pickle.load(file)
    file.close()
    
    ref_TaxID = file_dict['refTaxID']
    refA = rw_a2m.read_refseq(file_dict['fastaA'])
    refB = rw_a2m.read_refseq(file_dict['fastaB'])
    
    ######################        params        ###############################    
    LEVEL = ['genus','family','order','class','phylum','kingdom','Bacteria',
             'Eukaryota','Archaea','Viruses']
    DOMAIN = ['Bacteria','Eukaryota','Archaea','Viruses']

    
    ######################     1. read msa      ###############################
    msaA_data = rw_a2m.read_a2m(file_dict['msaA'],len(refA[-1][-1]),cov)
    print(f'1.1A sequence count: {len(msaA_data)}')

    msaA_data.pop(0)
    msaA_data = rw_a2m.parse_msa(msaA_data)
    print(f'1.1A init filter, sequence count: {len(msaA_data)}')

                   ############################################################
    msaB_data = rw_a2m.read_a2m(file_dict['msaB'],len(refB[-1][-1]),cov)
    print(f'1.1B sequence count: {len(msaB_data)}')

    msaB_data.pop(0)
    msaB_data = rw_a2m.parse_msa(msaB_data)
    print(f'1.1B init filter, sequence count: {len(msaB_data)}')


    ######################     2. common Tax     ##############################
    common_species,TaxID2common = cs.common_Tax(msaA_data,msaB_data,tax2id,tax_db)
    TaxID_dict = cs.Tax_groupmsa(common_species,TaxID2common,
                                 msaA_data,msaB_data,
                                 tax2id,tax_db)
    print(f'\n1.2 common tax count: {len(TaxID_dict)}')
    unpaired_num = cs.unpairedseq(TaxID_dict)
    print(f'1.2 unpaired count: {unpaired_num}')


    #####################     3. common rank     ##############################
    class_rank = cs.class_Tax(common_species,ref_TaxID,tax2id,tax_db)

    commomTax_msa = {}
    for subclass in LEVEL:
        commomTax_msa[subclass] = {}
        for TaxID in class_rank[subclass]:
            commomTax_msa[subclass][TaxID] = TaxID_dict[TaxID]

    print('\n1.3 common rank')
    for rank in class_rank:
        print(f'{rank:15s} {len(class_rank[rank])}')
        
            ###################################################################
    selectTax_A = []
    selectTax_B = []

    ref_domain = tax_db[ref_TaxID][-1][-1]
    
    if ref_domain in DOMAIN:
        DOMAIN.remove(ref_domain)
        select_level = LEVEL[:6]+[ref_domain]+DOMAIN
    else:
        select_level = LEVEL

    max_common_level = file_dict['max_common_level']
    if max_common_level == 'domain':
        select_level = select_level[:-3]
    elif max_common_level == 'life':
        pass
    else:
        temp = LEVEL.index(max_common_level)
        select_level = select_level[:(temp+1)]
    print(f'select level: {"  ".join(select_level)}')

    for subclass in select_level:
        for species in commomTax_msa[subclass].keys():
            selectTax_A += commomTax_msa[subclass][species][0]
            selectTax_B += commomTax_msa[subclass][species][1]

    msaA_data = selectTax_A
    msaB_data = selectTax_B

    common_species,TaxID2common = cs.common_Tax(msaA_data,msaB_data,tax2id,tax_db)
    TaxID_dict = cs.Tax_groupmsa(common_species,TaxID2common,
                                  msaA_data,msaB_data,
                                  tax2id,tax_db)
    print(f'1.3 common level filter,common tax count: {len(TaxID_dict)}')
    unpaired_num = cs.unpairedseq(TaxID_dict)
    print(f'1.3 unpaired count: {unpaired_num}')     
    TaxID_dict = cs.sorted_sim(TaxID_dict,refA,refB)

        
    ##################        4. final paired        ##########################     
    paired = final_pair(TaxID_dict,topn)
    paired_file = os.path.join(file_dict['outpath'],'paired.a3m')
    
    with open(paired_file,'w') as f:
        header = refA[-1][0].strip()+'&'+refB[-1][0].strip()
        seq = refA[-1][1].strip()+refB[-1][1].strip()
        f.write(header)
        f.write('\n')
        f.write(seq)
        f.write('\n')
            
        for header,seq in paired:
            f.write(header)
            f.write('\n')
            f.write(seq)
            f.write('\n')


faA,faB,msaA,msaB,output_dir,tax2id_file,taxdb_file,refTaxID,max_common_level,cov,topn = sys.argv[1:]

file_dict = {'fastaA':faA,
              'fastaB':faB,
              'msaA':msaA,
              'msaB':msaB,
              'outpath':output_dir,
              'tax2id':tax2id_file,
              'tax_db':taxdb_file,
              'refTaxID':int(refTaxID),  
              'max_common_level':max_common_level}

main(file_dict,float(cov)/100,int(topn))


















