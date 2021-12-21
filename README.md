# PPI-MSA_taxonomy-rank
This repositoriy provides the simplified approach described in XXX to generate the MSA of interalogs for intra-species protein-protein interactions, which also allows the user to restrict to different taxonomic ranks of the species of target PPI in the MSA generation.
The generated MSA of interologs can be further used as the input of Alphafold2 for protein complex structure prediction or as the input of DRN-1D2Dfor protein contact prediction.


## Installation
### 1.
    git clone https://github.com/ChengfeiYan/PPI-MSA_taxonomy-rank.git
### 2. Download the packaged species database
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
    8.  refTaxID: NCBI Taxonomy ID of target A&B (9506, ...).
    9.  max_common_level: selected rank (species, genus, family, order, class, phylum, kingdom, domain, life).
        “life” was used to represent all species in the UniRef100 database.
    10. cov: minimun converage with target sequence (50,60, ...).
    11. topn: In each species, the topn matching sequence is selected (1,2,3, ...).

### 2. Usage
    python pairing.py faA faB msaA msaB output_dir tax2id_file taxdb_file refTaxID max_common_level cov topn
