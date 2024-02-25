import re
import numpy as np
import pandas as pd
import pyslow5
import pysam
from tqdm import tqdm
from nanoCEM.normalization import normalize_signal,normalize_signal_with_lim
from nanoCEM.cem_utils import generate_bam_file,identify_file_path,prepare_move_table_file
# from nanoCEM.plot import draw_signal
import os
import argparse
score_dict={}
nucleotide_type=None
def extract_feature(line,strand,sig_move_offset,norm=True):
    pbar.update(1)
    read_id = line[0]
    # if read_id !='562eeb47-2b86-4fc7-abfc-5dce62f511ed':
    #     return None
    if read_id not in info_dict:
        return None
    # tackle moves tag
    moves_string = line[12]
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
    # assert len_raw_signal in paf and blow5

    if read is None :
        return None

    try:
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
        signal,shift,scale = normalize_signal_with_lim(signal)

    # index query and reference
    aligned_pair=info_dict[read_id]['pairs']
    qlen = info_dict[read_id]['query_length']
    try:
        assert qlen == len(event_length) + sig_move_offset
    except Exception:
        print("Warning: 1 read's length of event is not equal between bam and paf file")
        return None
    # correct index about DNA and RNA
    if (nucleotide_type == 'RNA' and strand == '+') or (nucleotide_type == 'DNA' and strand == '-'):
        aligned_pair[0] = qlen - aligned_pair[0] - 1


    if aligned_pair.shape[0]==0:
        return None
    read_pos = aligned_pair[0].values
    ref_pos = aligned_pair[1].values

    # extract raw signal by event length and event start
    total_feature_per_reads = []
    raw_signal_every = [signal[event_starts[x]:event_starts[x] + event_length[x]] for x in
                        read_pos]
    # calculate mean median and dwell time
    for i, element in enumerate(raw_signal_every):
        if len(element) == 0:
            continue
        temp = [read_id,np.mean(element), np.std(element), np.median(element), event_length[read_pos[i]],ref_pos[i]]
        total_feature_per_reads.append(temp)
    return total_feature_per_reads

def extract_pairs_pos(bam_file,position,length,chromosome,strand):

    result_dict={}
    for read in bam_file.fetch(chromosome,position-length, position+length+1):
        if strand == '+' and read.is_reverse:
            continue
        if strand == '-' and not read.is_reverse:
            continue
        # if read.qname == '7dcec2bc-65f3-41fe-8981-dd3af9fa6e67':
        #     print(1)
        start_position=read.reference_start
        end_position=read.reference_end
        if position < start_position or position > end_position:
            continue
        # unit
        aligned_pair = pd.DataFrame(np.array(read.aligned_pairs)).dropna().reset_index(drop=True)
        aligned_pair = aligned_pair[(aligned_pair[1] >= position - length) & (aligned_pair[1] <= position + length)]

        temp={}
        temp['pairs']=aligned_pair
        temp['query_length'] = read.query_length
        result_dict[read.qname] = temp
    return result_dict



def read_basecall_bam(path,position,reference,length,chrom,strand,sig_move_offset,kmer_length,subsample_ratio=1,norm=True,cpu=4,rna=True):
    global info_dict, s5,pbar,nucleotide_type
    slow5_file = path + ".blow5"
    bam_file = path + ".bam"
    identify_file_path(bam_file)
    identify_file_path(slow5_file)
    align_bam_file, paf_file = prepare_move_table_file(bam_file,reference,cpu,str(sig_move_offset),str(kmer_length))

    if rna:
        nucleotide_type ='RNA'
    else:
        nucleotide_type = 'DNA'
    bam_file = pysam.AlignmentFile(align_bam_file,'rb')

    info_dict=extract_pairs_pos(bam_file,position,length,chrom,strand)
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
    df["feature"] = df.apply(extract_feature,strand=strand,norm=norm,sig_move_offset=sig_move_offset,axis=1)
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
    # if rna_mode:
    #     if strand == '+':
    #         final_feature['position']=final_feature['position'] - (kmer_model-1)
    #     else:
    #         final_feature['position']=final_feature['position'] + (kmer_model-1)

    final_feature['Position'] = final_feature['Position'].astype(int).astype(str)
    print('Extracted ', num_aligned, ' aligned reads from blow5 files')

    # if num_aligned>50:
    #     dwell_filter_pctls = (0.5, 99.5)
    #     dwell_min, dwell_max = np.percentile(final_feature['Dwell time'].values, dwell_filter_pctls)
    #     final_feature = final_feature[(final_feature['Dwell time'] > dwell_min) & (final_feature['Dwell time'] < dwell_max)]

    return final_feature,num_aligned

