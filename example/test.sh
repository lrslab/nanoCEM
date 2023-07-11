## read tombo pipeline
cd data/ivt/
guppy_basecaller -i single/ -s ./guppy_out -c rna_r9.4.1_70bps_hac.cfg --recursive --device auto
cat guppy_out/*/*.fastq>final.fastq
tombo preprocess annotate_raw_with_fastqs --fast5-basedir single/ --fastq-filenames final.fastq --processes 16 --overwrite
tombo resquiggle single/ ../23S_rRNA.fasta --processes 16 --num-most-common-errors 5 --overwrite
# repeat preprocess on another sample
cd ../../
cd data/wt/
guppy_basecaller -i single/ -s ./guppy_out -c rna_r9.4.1_70bps_hac.cfg --recursive --device auto
cat guppy_out/*/*.fastq>final.fastq
tombo preprocess annotate_raw_with_fastqs --fast5-basedir single/ --fastq-filenames final.fastq --processes 16 --overwrite
tombo resquiggle single/ ../23S_rRNA.fasta --processes 16 --num-most-common-errors 5 --overwrite
cd ../../

current_events_magnifier.py tombo -i data/wt/single -c data/ivt/single -o tombo_result \
--chrom NR_103073.1 --strand + \
--pos 2030 \
--ref data/23S_rRNA.fasta \
--rna --cpu 4 --norm

## read f5c pipeline
cd data/ivt/
slow5tools f2s single -d blow5_dir
slow5tools merge blow5_dir -o file.blow5
slow5tools index file.blow5
minimap2 -ax map-ont -t 16 --MD ../23S_rRNA.fasta final.fastq | samtools view -hbS -F 260 - | samtools sort -@ 6 -o file.bam
samtools index file.bam
f5c resquiggle -c final.fastq file.blow5 -o file.paf
cd ../../
# repeat preprocess on another sample
cd data/wt/
slow5tools f2s single -d blow5_dir
slow5tools merge blow5_dir -o file.blow5
slow5tools index file.blow5
minimap2 -ax map-ont -t 16 --MD ../23S_rRNA.fasta final.fastq | samtools view -hbS -F 260 - | samtools sort -@ 6 -o file.bam
samtools index file.bam
f5c resquiggle -c final.fastq file.blow5 -o file.paf
cd ../../

current_events_magnifier.py f5c -i data/wt/file -c data/ivt/file -o f5c_result \
--chrom NR_103073.1 --strand + \
--pos 2030 \
--ref data/23S_rRNA.fasta \
--base_shift 2 --rna --norm
