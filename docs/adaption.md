# Adaption for preprocessing methods

## Introduction
nanoCEM now supports four alignment methods for nanopore signal: `f5c resquiggle`, `f5c eventalign`, `move_table`, and `tombo resquiggle`.

While `tombo resquiggle` and `f5c eventalign` align the signals to the reference sequence,
`f5c resquiggle` and `move_table` are alignment methods applied to the basecalled sequence.
Thus,these basecalled sequence methods require additional adaptions to index the reference sequence.

## Basecalled sequence methods adaption

### Signal index file acquisition
For `f5c resquiggle`, a **PAF** (Pairwise Alignment Format) file is generated, 
which records the alignment information and indices corresponding to the basecalled sequence. 

On the other hand, the `move_table` is stored in the **BAM** file generated from the basecalled sequence, 
and the related decoding and indexing is [more complex](https://github.com/hiruna72/squigualiser/blob/main/docs/move_table.md).
To address this issue, we utilized the `resquigualiser reform` tool available in the repository [**here**](https://github.com/hiruna72/squigualiser/blob/main/docs/reform.md)  to generate the `PAF` file.

### Re-indexing on the reference sequence
For such methods, we align the basecalled sequence to the reference by using the **CIGAR values** from the **BAM** file to generate an **index table**. 
This process involves discarding **insertions** and **deletions** while preserving **mismatches** to accurately map the basecalled sequence onto the reference.

<center>![adaption](adaption.bmp "adaption") </center>

## Base shift (only for f5c)

For `f5c` (both for `resquiggle and eventalign`), during nanopore signal alignment, it utilizes a **k-mer model** and assigns the alignment result to the **last base**. 
However, this is not the case for `tombo resquiggle` and `move_table`. 

To make their results as consistent as possible and enable comparisons, 
we introduced `--base_shift` option to align the result closer to the **middle** of the k-mer. 
For example, in R10 DNA sequence whose k-mer is 9-mer, if the basecalled sequence is ACACTACA**C** (9 nt), `f5c` will return only 1 event index of the last **C**.
But after turn on the **base_shift**, it will be shifted to the middle **T**(ACAC**T**ACAC)