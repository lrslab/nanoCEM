import os
import time
from collections import OrderedDict

import numpy as np
import pandas as pd


def read_fasta_to_dic(filename):
    """
    function used to parser small fasta
    still effective for genome level file
    """
    fa_dic = OrderedDict()

    with open(filename, "r") as f:
        for n, line in enumerate(f.readlines()):
            if line.startswith(">"):
                if n > 0:
                    fa_dic[short_name] = "".join(seq_l)  # store previous one

                full_name = line.strip().replace(">", "")
                short_name = full_name.split(" ")[0]
                seq_l = []
            else:  # collect the seq lines
                if len(line) > 8:  # min for fasta file is usually larger than 8
                    seq_line1 = line.strip()
                    seq_l.append(seq_line1)

        fa_dic[short_name] = "".join(seq_l)  # store the last one
    return fa_dic


def identify_file_path(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError("File do not exist! Please check your path : " + file_path)


def generate_bam_file(fastq_file, reference, cpu):
    bam_file = os.path.dirname(fastq_file) + '/' + os.path.basename(fastq_file).split('.')[0] + '.bam'
    if not os.path.exists(bam_file):
        cmds = 'minimap2 -ax map-ont -t ' + cpu + ' --MD ' + reference + ' ' + fastq_file + ' | samtools view -hbS -F ' + str(
            260) + '  - | samtools sort -@ ' + cpu + ' -o ' + bam_file
        print('Start to alignment ...')
        os.system(cmds)
        print('bam file is saved in ' + bam_file)
    else:
        print(bam_file + ' existed! Will skip the minimap2 ... ')
    cmds = 'samtools index ' + bam_file
    os.system(cmds)
    return bam_file


def run_samtools(fastq_file, location, reference, result_path, group, cpu):
    bam_file = generate_bam_file(fastq_file, reference, cpu)
    cmds = 'samtools mpileup ' + bam_file + ' -r ' + location + ' --no-output-ins --no-output-del -B -Q 0 -f ' + reference + ' -o ' + result_path + 'temp.txt'
    os.system(cmds)
    temp_file = pd.read_csv(result_path + 'temp.txt', sep='\t', header=None)
    temp_file['group'] = group
    os.system('rm ' + result_path + 'temp.txt')
    return temp_file


def build_out_path(results_path):
    if not os.path.exists(results_path):
        os.mkdir(results_path)
    else:
        print("Output file existed! It will be overwrite or add content after 5 secs ...")
        time.sleep(5)
        print("Continue ...")


def reverse_fasta(ref):
    base_dict = {'A': 'T', "T": "A", "C": "G", "G": "C"}
    ref = np.array(list(ref.upper()))
    for index, element in enumerate(ref):
        ref[index] = base_dict[element]
    reverse_fasta = list(reversed(''.join(ref)))
    return ''.join(reverse_fasta)


def save_fasta_dict(fasta_dict, path):
    f = open(path, 'w+')
    for key, value in fasta_dict.items():
        f.write('>' + key + '\n')
        line = len(value) // 80 + 1
        for i in range(0, line):
            f.write(value[i * 80:(i + 1) * 80] + '\n')
        f.close()
# fasta=read_fasta_to_dic("../example/data/23S_rRNA.fasta")
# for key,value in fasta.items():
#     fasta[key]=reverse_fasta(value)
# save_fasta_dict(fasta,'../example/data/23S_rRNA_re.fasta')
