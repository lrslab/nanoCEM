from collections import OrderedDict
import numpy as np

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

def reverse_fasta(ref):
    base_dict={'A':'T',"T":"A","C":"G","G":"C"}
    ref=np.array(list(ref.upper()))
    for index,element in enumerate(ref):
        ref[index] = base_dict[element]
    reverse_fasta = list(reversed(''.join(ref)))
    return ''.join(reverse_fasta)

def save_fasta_dict(fasta_dict,path):
    f=open(path,'w+')
    for key,value in fasta_dict.items():
        f.write('>'+key+'\n')
        line = len(value) //80 + 1
        for i in range(0, line):
            f.write(value[i*80:(i+1)*80]+'\n')
        f.close()
# fasta=read_fasta_to_dic("/data/Ecoli_23s/23S_rRNA.fasta")
# for key,value in fasta.items():
#     fasta[key]=reverse_fasta(value)
# save_fasta_dict(fasta,'test.fasta')