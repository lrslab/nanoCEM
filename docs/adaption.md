# ![logo](logo_tiny.png "nanoCEM") Adaption for preprocessing methods

## Introduction
nanoCEM now supports four alignment methods for nanopore signal: `f5c resquiggle`, `f5c eventalign`, `move_table` from basecaller, and `tombo resquiggle`.

While `tombo resquiggle` and `f5c eventalign` align the signals to the reference sequence,
`f5c resquiggle` and `move_table` from basecaller are alignment methods applied to the basecalled sequence.
Thus,these basecalled sequence methods require additional adaptions to index the reference sequence.

## Basecalled sequence methods adaption

### Signal index file acquisition
For `f5c resquiggle`, a **paf** (Pairwise Alignment Format) file is generated, 
which records the alignment information and indices corresponding to the basecalled sequence. 

On the other hand, the `move_table` is stored in the **bam** file generated from the basecaller, 
for the related decoding and indexing methods, detailed description is [here](https://github.com/hiruna72/squigualiser/blob/main/docs/move_table.md).
To address this issue, we utilized the `squigualiser reform` tool available in the repository [**here**](https://github.com/hiruna72/squigualiser/blob/main/docs/reform.md)  to generate the `PAF` file.

### Re-indexing on the reference sequence
For such methods, we align the basecalled sequence to the reference by using the **CIGAR values** from the **bam** file to generate an **index table**. 
This process involves discarding **insertions** and **deletions** while preserving **mismatches** to accurately map the basecalled sequence onto the reference.

<center>![adaption](adaption.bmp "adaption") </center>

## Base shift (only for f5c)

For `f5c` (both for `resquiggle and eventalign`), during nanopore signal alignment, it utilizes a **k-mer model** and assigns the alignment result to the **first base**. However, this is not the case for `tombo resquiggle` and `move_table`.
Their developers noticed this problem and applied a complicated strategy to tackle it in this [link](https://github.com/hiruna72/squigualiser/blob/main/docs/base_shift_and_eventalignment.md).

But in nanoCEM, to make their results as consistent as possible and enable comparisons with them, 
we applied a simple strategy to shift f5c's result and introduced `--base_shift` option to align the result closer to the **middle** of the k-mer. 

For example, in R10 DNA sequence whose k-mer is 9-mer, if the basecalled sequence is **A**CACTACAC (9 nt), `f5c` will return only 1 event index of the first **A**.
But after turn on the **base_shift** mode, it will be shifted to the middle **T**(ACAC**T**ACAC)


| kit    | type | k-mer model  |base shift number |example |
|--------|----------|-------|-------|-------|
| R9    | DNA | 6 |2 |**C**TACAC → CT**A**CAC |
| R10    | DNA | 9 |4 |**A**CACTACAC → ACAC**T**ACAC |
| R9    | RNA002 | 5 |2 |**U**ACAC → UA**C**AC |
| R9    | RNA004 | 9 |4 |**A**CACUACAC → ACAC**U**ACAC |
