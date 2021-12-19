#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 12:06:13 2021

@author: yunda_si
"""

import rw_msa
import numpy as np
import copy



def cal_similarity(square_mat,np_fasta):

    for i in range(len(square_mat)):
        idendity = np.sum(np_fasta[i]==np_fasta[-1])
        square_mat[i] = idendity
        
    return square_mat


def Tax_groupmsa(common_species,TaxID2common,
                 parsed_msaA,parsed_msaB,
                 tax2id,tax_db):

    TaxID_dict = {}
    for TaxID in common_species:
        TaxID_dict[TaxID] = [[],[]]

    for fasta in parsed_msaA:

        TaxID = fasta[0][4]
        if TaxID not in TaxID2common:
            continue

        TaxID_rank = tax_db[TaxID][0]
        if TaxID_rank == 'species':
            TaxID_dict[TaxID][0].append(fasta)
        else:
            speciesid = tax2id[tax_db[TaxID][-1][2]]
            TaxID_dict[speciesid][0].append(fasta)

    for fasta in parsed_msaB:

        TaxID = fasta[0][4]
        if TaxID not in TaxID2common:
            continue

        TaxID_rank = tax_db[TaxID][0]
        if TaxID_rank == 'species':
            TaxID_dict[TaxID][1].append(fasta)
        else:
            speciesid = tax2id[tax_db[TaxID][-1][2]]
            TaxID_dict[speciesid][1].append(fasta)

    return TaxID_dict



def common_Tax(parsed_msaA,parsed_msaB,tax2id,tax_db):

    s2s_dict = {}  #species2strain_dict
    ####################            common Tax            #####################
    TaxID_msaA = set()
    for fasta in parsed_msaA:
        TaxID = fasta[0][4]
        try:
            rank_self = tax_db[TaxID][0]
            if rank_self == 'species':
                if set(tax_db[TaxID][-1][3:]) != {''}:
                    TaxID_msaA.add(TaxID)

            elif rank_self == 'strain':
                species_TaxID = tax2id[tax_db[TaxID][-1][2]]
                if set(tax_db[species_TaxID][-1][3:]) != {''}:
                    TaxID_msaA.add(species_TaxID)

                    if species_TaxID not in s2s_dict.keys():
                        s2s_dict[species_TaxID] = []
                    s2s_dict[species_TaxID].append(TaxID)
            else:
                pass
        except:
            pass
    TaxID_msaA.discard(-1)

    TaxID_msaB = set()
    for fasta in parsed_msaB:
        TaxID = fasta[0][4]
        try:
            rank_self = tax_db[TaxID][0]
            if rank_self == 'species':
                if set(tax_db[TaxID][-1][3:]) != {''}:
                    TaxID_msaB.add(TaxID)

            elif rank_self == 'strain':
                species_TaxID = tax2id[tax_db[TaxID][-1][2]]
                if set(tax_db[species_TaxID][-1][3:]) != {''}:
                    TaxID_msaB.add(species_TaxID)

                    if species_TaxID not in s2s_dict.keys():
                        s2s_dict[species_TaxID] = []
                    s2s_dict[species_TaxID].append(TaxID)
            else:
                pass
        except:
            pass
    TaxID_msaB.discard(-1)

    common_species = TaxID_msaA & TaxID_msaB

    ####################            add strain            #####################
    s2s_set = set(s2s_dict.keys())
    common = common_species & s2s_set

    TaxID2common = copy.deepcopy(common_species)
    if common:
        for species in common:
            TaxID2common = TaxID2common|set(s2s_dict[species])

    return common_species,TaxID2common



def class_Tax(common_TaxID,ref_TaxID,tax2id,tax_db):

    level = ['genus','family','order','class','phylum','kingdom','Archaea',
             'Bacteria','Eukaryota','Viruses']
    
    class_rank = {rank:[] for rank in level}

    reftax_rank = list(tax_db[ref_TaxID][-1][3:])
    reftax_rank = ['refgap' if i == '' else i for i in reftax_rank]
    refid_rank = [tax2id.get(i,-3) for i in reftax_rank]

    for TaxID in common_TaxID:

        Tax_rank = list(tax_db[TaxID][-1][3:])
        TaxID_rank = [tax2id.get(i,-2) for i in Tax_rank]

        for index,rank in enumerate(level[:6]):
            have_common = TaxID_rank[index] == refid_rank[index]
            if have_common:
                class_rank[rank].append(TaxID)
                break

        if index == 5 and not have_common:
            if Tax_rank[index+1] ==  'Archaea':
                class_rank['Archaea'].append(TaxID)

            elif Tax_rank[index+1] ==  'Bacteria':
                class_rank['Bacteria'].append(TaxID)

            elif Tax_rank[index+1] ==  'Eukaryota':
                class_rank['Eukaryota'].append(TaxID)

            elif Tax_rank[index+1] ==  'Viruses':
                class_rank['Viruses'].append(TaxID)             
            else:
                print(f'{TaxID} no class information')
                
    return class_rank



def unpairedseq(TaxID_dict):

    count = 0
    for TaxID in TaxID_dict:
        faA_list,faB_list = TaxID_dict[TaxID]
        faA_num = len(faA_list)
        faB_num = len(faB_list)
        
        if faA_num+faB_num>2:
            count += min(faA_num,faB_num)
    
    return count


def sorted_sim(TaxID_dict,seqA,seqB):

    for TaxID in TaxID_dict:
        faA_list,faB_list = TaxID_dict[TaxID]

        faA_num = len(faA_list)
        faB_num = len(faB_list)

        if faA_num>1:
            faA_list.append(seqA)
            np_fastaA = rw_a2m.encode_a2m(faA_list)
            square_mat = np.zeros((faA_num),dtype=np.float64)
            square_mat = cal_similarity(square_mat,np_fastaA)
            filtered_index = square_mat.argsort()[::-1]

            sorted_faA_list = []
            for index in filtered_index:
                sorted_faA_list.append(faA_list[index])
            TaxID_dict[TaxID][0] = sorted_faA_list

        if faB_num>1:
            faB_list.append(seqB)
            np_fastaB = rw_a2m.encode_a2m(faB_list)
            square_mat = np.zeros((faB_num),dtype=np.float64)
            square_mat = cal_similarity(square_mat,np_fastaB)
            filtered_index = square_mat.argsort()[::-1]

            sorted_faB_list = []
            for index in filtered_index:
                sorted_faB_list.append(faB_list[index])
            TaxID_dict[TaxID][1] = sorted_faB_list

    return TaxID_dict

