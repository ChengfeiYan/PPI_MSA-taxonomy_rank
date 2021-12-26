### Some changes to the alphafold2 code:
Change "predicted_aligned_error" to "predicted_aligned_error_2d" in line 106 of alphafold/common/confidence.py ", 
then the "predicted_aligned_error" in the "result.pkl" file output by AF2 corresponds to the original data used to 
calculate the ptm, "predicted_aligned_error_2d" corresponds to the original PAE matrix.

Then input "result.pkl" into ptm.py and calculate the ptm value after removing the linker.

### Parameters of the ptm.py script：
1. af2_path： Installation path for alphafold2.
2. result.pkl: Location of alphafold2 prediction result.
3. lenA: Length of sequence A.
4. len_gap: Length of gap.

### Usage：
    python ptm.py af2_path result_pkl lenA len_gap
