# ![logo](logo_tiny.png "nanoCEM")Command line arguments
We provide 3  scripts to run the analysis and visualize nanopore data as the following.

alignment_magnifier
=============


| Argument name                       | Required | Description                                     |
|---------------------------------|----------|-------------------------------------------------|
| -i, --input_fastq               | Yes       | input FASTQ file                                |
| -c, --control_fastq             | No       | control FASTQ file, will run as single-mode if missing            |
| -o, --output                    | No       | output folder's name (Default:nanoCEM_result)    |
| --chrom                         | Yes       | gene/transcript or chromosome name (header of your reference file)|
| --pos                           | Yes       | Site of your interest (**1-based**)         |
| -r, --ref                       | Yes       | reference file (fasta/fna/fa)          |
| --len                           | No       | region around the position   (Default:10)|
| -t,--cpu                           | No       | number of processes   (Default:4)          |
| --strand                        | No       | strand of your interest  (Default:+)          |
| --rna                           | No       | turn on the RNA mode (Default:False)       |
| -h, --help                      | No       | show this help message and exit          |

current_events_magnifier 
=============
These are some differences in parameters when using "current_events_magnifier" with `tombo` and `f5c` tools.
For tombo, you can view it for 

    current_events_magnifier tombo -h

|  Argument  name            | Required | Description                                                     |
|---------------------------------|----------|-----------------------------------------------------------------|
| -i , --fast5          | Yes       | the input single-format fast5's folder                            |
| -c , --control_fast5 | No       | the control single-format fast5's folder, will run as single-mode if missing              |
| -o , --output       | No       | output_file (Default:nanoCEM_result)                                                     |
| --ref                        | Yes       |reference file (fasta/fna/fa)          |
| --chrom                    | Yes       | gene/transcript or chromosome name (headdf of your fasta file)               |
| --pos                        | Yes       | site of your interest        (**1-based**)          |
| --len                        | No       | region around the position (Default:10)                          |
| -t, --cpu                | No       | num of process (Default:8)                                       |
| --subsample-ratio | No       |  ratio of reads will be used to plot (Default:1)                 |
| --strand                  | No       | strand of your interest (Default:+)   |
| --rna                           | No       | turn on the RNA mode       (Default:False)                                      |
| --norm                          | No       | turn on the normalization      (Default:False)                  |
| --basecall_group                | No       | re-squiggle result index information into.(Default:RawGenomeCorrected_000)           |
| --basecall_subgroup             | No       | basecall sequence inforamtion into.  (Default:BaseCalled_template) |
| -h, --help                      | No       | Show this help message and exit                                  |

And f5c
    current_events_magnifier f5c -h

|  Argument name            | Required | Description                                                     |
|---------------------------------|----------|-----------------------------------------------------------------|
| -i , --input          | Yes       | the input files' prefix                            |
| -c , --control | No       | the control files' prefix, will run as single-mode if missing           |
| -o , --output       | No       | output_file (Default:nanoCEM_result)               |
| --ref                        | Yes       |reference file (fasta/fna/fa)          |
| --chrom                    | Yes       | gene/transcript or chromosome name (header of your reference file)               |
| --pos                        | Yes       | site of your interest   (**1-based**)          |
| --len                        | No       | region around the position (Default:10)                          |
| -t, --cpu                | No       | num of process (Default:8)                                       |
| --subsample-ratio | No       |  ratio of reads will be used to plot (Default:1)                 |
| --strand                  | No       | strand of your interest (Default:+)   |
| --rna                           | No       | turn on the RNA mode       (Default:False)                                      |
| --norm                          | No       | turn on the normalization      (Default:False)                  |
| --base_shift                          | No       |  for shifting the result from f5c      (Default:2)                  |
| --pore | No       |  type of pore choose from `r9` and `r10` (Default: `r9`)                 |
| -h, --help                      | No       | Show this help message and exit                                  |

extract_sub_fast5_from_bam
==================
 To extract a subset of a single-format fast5 and saved in new folder.
 If you do not input an interest region, all the reads that have been aligned will be saved in the new folder.

| Argument name| Required | Description                                       |
|-------------------|----------|---------------------------------------------------|
| -f , --fast5  | Yes  | path to single-format fast5 files                           |
| -b , --bam  | Yes    | path to bam file                                  |
| -o , --output  | No  | output file   (Default:subsample_single)        |
| -t CPU, --cpu  | No       | number of processes  (Default:4)    |
| --chrom      | No       | gene/transcript or chromosome name (header of your reference file)|
| --pos          | No       | site of your interest    (**1-based**)         |
| -h, --help        | No       | show this help message and exit                    |