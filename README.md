# PPI-MSA_taxonomy-rank
This repositoriy provides the code implementation of the algorithm described in [1] to generate the MSA of interalogs for intra-species protein-protein interactions, which also allows the user to restrict interologs to different taxonomic ranks of the species of target PPI in the MSA generation.
The generated MSA of interologs can be further used as the input of [Alphafold2](https://github.com/deepmind/alphafold) for protein complex structure prediction or as the input of [DRN-1D2D](https://github.com/ChengfeiYan/DRN-1D2D) for protein contact prediction.

In addition, you can use /ptm/ptm.py to calculate the ptm value after removing the linker.

## Installation
### 1.
    git clone https://github.com/ChengfeiYan/PPI-MSA_taxonomy-rank.git
### 2. Download the packaged species database (tax2id.pkl and taxdict.pkl)
    https://drive.google.com/file/d/1omdyzewWyx7E-orK6yu6wtMskfMPUmsj/view?usp=sharing

## Usage
### 1. params
    1.  faA: fasta file corresponding to target A.
    2.  faB: fasta file corresponding to target B.
    3.  msaA: a3m file corresponding to target A.
    4.  msaB: a3m file corresponding to target B.
    5.  output_dir: /
    6.  tax2id_file: Python's dictionary file, key is the name of the species, value is the NCBI Taxonomy ID  
            corresponding to the name.
    7.  taxdict_file: Python's dictionary file, key is NCBI Taxonomy ID, value is the taxonomic ranks
            (species, genus, family, order, class, phylum, kingdom, domain). 
    8.  refTaxID: NCBI Taxonomy ID of target A&B,Once you know the name of the species, 
            you can get the NCBI Taxonomy ID corresponding to the name from https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi.
    9.  max_common_level: selected rank (species, genus, family, order, class, phylum, kingdom, domain, life).
            “life” was used to represent all species in the UniRef100 database.
    10. cov: minimun converage with target sequence, the recommended value is 60.
    11. topn: In each species, the topn matching sequence is selected，When set to 1, it means that only the top pair is used, 
            and when set to 100000, it means that all pairs are used..

### 2. Usage
    python pairing.py faA faB msaA msaB output_dir tax2id_file taxdb_file refTaxID max_common_level cov topn

## Citing:
[1]. Protein Complex Structure Prediction Powered by Multiple Sequence Alignment of Interologs from Multiple Taxonomic Ranks and AlphaFold2
Yunda Si, Chengfei Yan. Briefings in Bioinformatics 23, bbac208, 2022.
If you meet any problem in installing or running the program, please contact chengfeiyan@hust.edu.cn.
