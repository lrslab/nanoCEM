import re
import numpy as np
import pandas as pd
import pyslow5
import pysam
from tqdm import tqdm
from nanoCEM.normalization import normalize_signal,normalize_signal_with_lim
from nanoCEM.cem_utils import generate_bam_file,identify_file_path,generate_paf_file_eventalign,base_shift_dict
# from nanoCEM.plot import draw_signal
# import os
# import argparse
score_dict={}
nucleotide_type=None

def extract_feature(line,strand,position,windows_length,base_shift=2,norm=True):
    global nucleotide_type
    pbar.update(1)
    read_id = line[0]
    # if read_id == '0a017d09-5c09-456d-8875-635f4c2c1380':
    #     print(1)
    if line[4] != strand:
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
            continue
        else:
            event_length.append(int(item))
    # build event_length from move table
    read = s5.get_read(read_id, aux=["read_number", "start_mux"],pA=True)
    start_index = line[2]
    end_index = line[3]
    event_length = np.array(event_length)

    # identify RNA or DNA
    if nucleotide_type is None:
        if line[7] > line[8]:
            nucleotide_type = 'RNA'
        else:
            nucleotide_type = 'DNA'
    #  assert len_raw_signal in paf and blow5
    try:
        # assert end_index-start_index == np.sum(event_length)
        assert np.abs(line[8] - line[7]) == len(event_length)
        assert read['len_raw_signal'] == line[1]
    except Exception:
        print("Warning: 1 read's length of signal is not equal between blow5 and paf")
        return None

    signal = read['signal']

    signal = signal[start_index:end_index]
    if norm:
        signal, shift, scale = normalize_signal_with_lim(signal)
    event_starts = event_length.cumsum()
    event_starts = np.insert(event_starts, 0, 0)[:-1]

    # base shift
    ref_start = np.min([line[7],line[8]])
    ref_end = np.max([line[7], line[8]])
    ref_start = ref_start - base_shift
    ref_end = ref_end - base_shift
    start_position = np.max([ref_start, position - windows_length]) - ref_start
    end_position = np.min([ref_end, position + windows_length]) - ref_start

    # index query and reference map index
    if (nucleotide_type == 'RNA' and strand=='+') or (nucleotide_type == 'DNA' and strand=='-'):
        end_pos = line[10] - start_position - 1
        start_pos = line[10] - end_position - 1
    else:
        end_pos = end_position
        start_pos = start_position

    end_pos = np.min([end_pos, line[10] - 1])
    # extract raw signal by event length and event start
    total_feature_per_reads = []
    raw_signal_every = [signal[event_starts[x]:event_starts[x] + event_length[x]] for x in
                        range(start_pos,end_pos+1)]
    if (nucleotide_type == 'RNA' and strand=='+') or (nucleotide_type == 'DNA' and strand=='-'):
        raw_signal_every.reverse()
    # calculate mean median and dwell time
    for i, element in enumerate(raw_signal_every):
        if len(element) == 0:
            continue
        temp = [read_id,np.mean(element), np.std(element), np.median(element), len(element), str(start_position+ref_start+i)]
        total_feature_per_reads.append(temp)
    return total_feature_per_reads

def extract_pairs_pos(bam_file,position,length,chromosome,strand):

    result_dict={}
    for read in bam_file.fetch(chromosome,position-length, position+length+1):
        if read.is_supplementary or read.is_secondary:
            continue
        if strand == '+' and read.is_reverse:
            continue
        if strand == '-' and not read.is_reverse:
            continue
        # if read.qname == 'db71b047-e073-42cf-833b-a3ccdd9459b3':
        #     print(1)
        start_position=read.reference_start
        end_position=read.reference_end
        if position < start_position or position > end_position:
            continue
        # unit
        temp={}
        temp['ref_length'] = read.reference_length
        result_dict[read.qname] = temp
    return result_dict

def read_blow5(path,position,reference,length,chrom,strand,pore,aligned_file,subsample_ratio=1,base_shift='auto',norm=True,cpu=4,rna=True):
    global s5,pbar
    slow5_file = path + ".blow5"
    fastq_file = path + ".fastq"
    identify_file_path(slow5_file)
    if rna:
        nucleotide_type = "RNA"
    else:
        nucleotide_type = 'DNA'
    if base_shift == 'auto':
        base_shift = base_shift_dict[pore+nucleotide_type+strand]
    else:
        base_shift = base_shift
    if aligned_file is None:
        identify_file_path(fastq_file)
        fastq_file, bam_file = generate_bam_file(fastq_file, reference, cpu, subsample_ratio)
        paf_file = generate_paf_file_eventalign(fastq_file,slow5_file,bam_file,reference,pore,rna,cpu)
    else:
        print("Detected that a PAF file already exists,will skip the eventalign step")
        bam_file = aligned_file+'.bam'
        paf_file = aligned_file+'.paf'
        identify_file_path(bam_file)
        identify_file_path(paf_file)
    bam_file = pysam.AlignmentFile(bam_file,'rb')

    info_dict = extract_pairs_pos(bam_file,position,length,chrom,strand)
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
    df["feature"] = df.apply(extract_feature,base_shift=base_shift,strand=strand,position=position,windows_length=length,norm=norm,axis=1)
    pbar.close()

    df.dropna(inplace=True)
    num_aligned = df.shape[0]
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

    return final_feature,num_aligned,nucleotide_type

