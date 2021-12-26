#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 20:26:21 2021

@author: yunda_si
"""

import sys
import os
import pickle
import numpy as np

af2_path, result_pkl, lenA, len_gap = sys.argv[1:]

sys.path.append(os.path.join(af2_path,'common'))
from  confidence import predicted_tm_score


result = pickle.load(open(result_pkl,'rb'))
pae_logits = result['predicted_aligned_error']['logits']
pae_breaks = result['predicted_aligned_error']['breaks']

lenA = int(lenA)
len_gap = int(len_gap) 

new_logits = np.delete(pae_logits,np.arange(lenA,lenA+len_gap),axis=0)
new_logits = np.delete(new_logits,np.arange(lenA,lenA+len_gap),axis=1)
ptm = predicted_tm_score(new_logits,pae_breaks)    

print(f'ptm: {ptm}')
   