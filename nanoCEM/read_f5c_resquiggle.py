import re
import numpy as np
import pandas as pd
import pyslow5
import pysam
from tqdm import tqdm
from nanoCEM.normalization import normalize_signal,normalize_signal_with_lim
from nanoCEM.cem_utils import generate_bam_file,identify_file_path,generate_paf_file_resquiggle,ker_model_size,caculate_base_shift_size
# from nanoCEM.plot import draw_signal
from read_move_table import  extract_pairs_pos
# import os
# import argparse
score_dict={}
nucleotide_type=None
def extract_feature(line,strand,kmer_model,base_shift,norm=True):
    global nucleotide_type
    pbar.update(1)
    read_id = line[0]
    # if read_id =='a3f64d91-bdd6-4744-8afc-eb9a9a1b00d6':
    #     print(1)
    if read_id not in info_dict:
        return None
    # tackle moves tag
    moves_string = line[14]
    moves_string = re.sub('ss:Z:', '', moves_string)
    moves_string = re.sub('D', 'D,', moves_string)
    moves_string = re.sub('I', 'I,', moves_string)
    # print(moves_string)
    moves = re.split(r',+', moves_string)
    moves = moves[:-1]
    # extract index and generate event_length and event_start
    insertion = 0
    event_length = []
    for i,item in enumerate(moves):
        if 'D' in item:
            deletion = int(item[:-1])
            for i in range(deletion):
                event_length.append(0)
        elif 'I' in item:
            if i == 0 :
                continue
            else:
                return None
        elif '=' in item:
            return None
        else:
            event_length.append(int(item))
    # build event_length from move table
    read = s5.get_read(read_id, aux=["read_number", "start_mux"],pA=True)
    start_index = line[2]
    end_index = line[3]
    event_length = np.array(event_length)

    # identify RNA or DNA
    if nucleotide_type is None:
        if line[7]>line[8]:
            nucleotide_type='RNA'
        else:
            nucleotide_type='DNA'
    #  assert len_raw_signal in paf and blow5
    try:
        assert end_index-start_index == np.sum(event_length)
        assert read['len_raw_signal'] == line[1]

    except Exception:
        print("Warning: 1 read's length of signal is not equal between blow5 and paf")
        return None

    signal = read['signal']

    # create event start and flip all table
    signal = signal[start_index:end_index]
    event_starts = event_length.cumsum()
    event_starts = np.insert(event_starts, 0, 0)[:-1]
    if norm:
        signal,_,_ = normalize_signal_with_lim(signal)

    # index query and reference
    aligned_pair=info_dict[read_id]['pairs']
    qlen = info_dict[read_id]['query_length']
    # assert len(event_length) + kmer_model - 1 == qlen

    # correct index about DNA and RNA
    if (nucleotide_type == 'RNA' and strand == '+') or (nucleotide_type == 'DNA' and strand == '-'):
        aligned_pair[0] = qlen - aligned_pair[0] - 1
        aligned_pair[0] = aligned_pair[0] - kmer_model + 1

    if base_shift != 0:
        aligned_pair[1] = aligned_pair[1] + base_shift


    if aligned_pair.shape[0]==0:
        return None
    read_pos = aligned_pair[0].values
    ref_pos = aligned_pair[1].values

    # extract raw signal by event length and event start
    total_feature_per_reads = []
    try:
        raw_signal_every = [signal[event_starts[x]:event_starts[x] + event_length[x]] for x in
                            read_pos]
        # if aligned_pair.shape[0] < 11:
        #     draw_signal(signal[event_starts[read_pos[0]]:event_starts[read_pos[-1]+1]], event_starts[read_pos[0]:read_pos[-1]+1]-event_starts[read_pos[0]], base_list)
    except Exception as e:
        print(1)
    # calculate mean median and dwell time
    for i, element in enumerate(raw_signal_every):
        if event_length[read_pos[i]] == 0:
            continue
        temp = [read_id,np.mean(element), np.std(element), np.median(element), event_length[read_pos[i]],ref_pos[i]]
        total_feature_per_reads.append(temp)
    return total_feature_per_reads


def read_blow5(path,position,reference,length,chrom,strand,pore,subsample_ratio=1,base_shift=True,norm=True,cpu=4,rna=True):
    global s5,pbar,info_dict
    slow5_file = path + ".blow5"
    fastq_file = path + ".fastq"
    identify_file_path(fastq_file)
    identify_file_path(slow5_file)

    if rna:
        nucleotide_type='RNA'
    else:
        nucleotide_type='DNA'
    kmer_model = ker_model_size[pore+'+'+nucleotide_type]
    if base_shift:
        base_shift = caculate_base_shift_size(kmer_model,strand)
    else:
        base_shift = 0

    fastq_file, bam_file = generate_bam_file(fastq_file, reference, cpu, subsample_ratio)
    paf_file = generate_paf_file_resquiggle(fastq_file,slow5_file,pore,rna,cpu)
    bam_file = pysam.AlignmentFile(bam_file,'rb')
    info_dict = extract_pairs_pos(bam_file,position,length+base_shift,chrom,strand)
    if info_dict == {}:
        raise RuntimeError("There is no read aligned on this position")
    info_df = pd.DataFrame(list(info_dict.keys()))

    s5 = pyslow5.Open(slow5_file, 'r')
    df=pd.read_csv(paf_file,sep='\t',header=None)
    df=pd.merge(df,info_df,how='inner',on=0)
    if df.shape[0] == 0:
        raise RuntimeError("cannot found the record from bam in your paf file. Please check your f5c command ... ")
    if df.shape[0] / info_df.shape[0] < 0.8:
        print('There are '+str(info_df.shape[0]-df.shape[0])+" reads not found in your paf file ...")
    pbar = tqdm(total=df.shape[0], position=0, leave=True)
    df["feature"] = df.apply(extract_feature,kmer_model=kmer_model,base_shift=base_shift,strand=strand,norm=norm,axis=1)
    pbar.close()

    df.dropna(inplace=True)
    num_aligned = df.shape[0]
    if subsample_ratio<1:
        df=df.sample(frac=subsample_ratio)
    final_feature=[]
    for item in df["feature"]:
        final_feature.extend(item)
    final_feature=pd.DataFrame(final_feature)
    final_feature.columns=['Read ID','Mean','STD','Median','Dwell time','Position']
    final_feature = final_feature[(final_feature['Position']>=position-length)&(final_feature['Position']<=position+length)]

    final_feature['Position'] = final_feature['Position'].astype(int).astype(str)
    print('Extracted ', num_aligned, ' aligned reads from blow5 files')

    # if num_aligned>50:
    #     dwell_filter_pctls = (0.5, 99.5)
    #     dwell_min, dwell_max = np.percentile(final_feature['Dwell time'].values, dwell_filter_pctls)
    #     final_feature = final_feature[(final_feature['Dwell time'] > dwell_min) & (final_feature['Dwell time'] < dwell_max)]

    return final_feature,num_aligned,nucleotide_type

