# Current_Magnifier
`Current_Magnifier` is a sample tool designed to visualize the features(Mean,Median,Dwell_time,STD) that distinguish between two groups of ONT data from the same species at the site level.
If you wanna view single read signal or raw signal, [Squigualiser](https://github.com/hiruna72/squigualiser) is recommended
## Still developing
## Quick start
### Preprocessing
If you used R10 or wanna use the last resquiggle program(f5c v1.2.0),
```sh
# if fast5 is single format 
# single is fast5s-base-directory
single_to_multi_fast5 -i single/ -s multi -n 1000 --recursive

slow5tools f2s multi/ -d blow5_dir
slow5tools merge blow5_dir -o file.blow5
slow5tools index file.blow5

minimap2 -ax map-ont -t 16 <reference-fasta> <reads.fastq> | samtools view -hbS -F 260 - | samtools sort -@ 6 -o file.bam
samtools index file.bam

f5c resquiggle -c final.fastq file.blow5 -o file.paf

python read_f5c_resquiggle.py -i file -c control -o f5c_result --chrom NC_000xxx --strand + --pos 3929 --len 5

```
If you used Tombo(v1.5.0),
```sh
# if fast5 is not single format 
# single is fast5s-base-directory
multi_to_single_fast5 -i single/ -s multi -n 1000 --recursive

tombo preprocess annotate_raw_with_fastqs --fast5-basedir  single/ --fastq-filenames <reads.fastq>

tombo resquiggle single. <reference-fasta> --processes 4 --num-most-common-errors 5

python read_tombo_resquiggle.py -i wt/single -c control/single -o tombo_result --chrom NC_000xxx --strand + --pos 3929 --len 5 --=cpu 4

```

