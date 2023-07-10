import numpy as np
import h5py
import pandas as pd
import plotnine as p9
import argparse
import os
import multiprocessing
from tqdm import tqdm
import random
# from ._stats import c_new_mean_stds
from nanoCEM.normalization import normalize_signal,normalize_signal_with_lim



def matrix_append(total_matrix,new_matrix,index):
    if new_matrix.shape[0]==0:
        return total_matrix
    len_new_matrix=new_matrix.shape[0]
    if index+len_new_matrix>total_matrix.shape[0]:
        temp_matrix=np.zeros((500000,5),dtype=float)
        total_matrix=np.vstack((total_matrix,temp_matrix))
    total_matrix[index:index+len_new_matrix]=new_matrix
    return total_matrix

def get_label_raw(fast5_fn, basecall_group, basecall_subgroup ,chromosome, position,strand):
    ##Open file
    # fast5_fn='/media/zhguo/QUO/ivt_positive/ivt_positive/single/42/878b9de7-043b-4d5e-b194-7771c5767ff5.fast5'
    # if '4a4438ea-fc46-42aa-9788-8d54f37471b9' not in fast5_fn:
    #     return None
    try:
        fast5_data = h5py.File(fast5_fn, 'r')
    except IOError:
        raise IOError('Error opening file. Likely a corrupted file.')

    # Get raw data
    try:
        raw_dat = list(fast5_data['/Raw/Reads/'].values())[0]
        raw_dat = raw_dat['Signal'][()]
    except:
        raise RuntimeError(
            'Raw data is not stored in Raw/Reads/Read_[read#] so ' +
            'new segments cannot be identified.')

    # Read corrected data
    try:
        corr_data = fast5_data['/Analyses/' + basecall_group + '/' + basecall_subgroup + '/Events']
        align_data = fast5_data['/Analyses/' + basecall_group + '/' + basecall_subgroup + '/Alignment']
        corr_attrs = dict(list(corr_data.attrs.items()))
        align_attrs=dict(list(align_data.attrs.items()))
        corr_data = corr_data[()]
    # ~ .value
    except:
        raise RuntimeError(('Corrected data not found.'))

    # fast5_info = fast5_data['UniqueGlobalKey/channel_id'].attrs
    # sampling_rate = fast5_info['sampling_rate'].astype('int_')

    # Reading extra information
    start_position=align_attrs['mapped_start']
    corr_start_rel_to_raw = corr_attrs['read_start_rel_to_raw']
    chrome=align_attrs["mapped_chrom"]
    if chrome != chromosome:
        return None
    if strand != align_attrs["mapped_strand"]:
        return None
    if position <= align_attrs['mapped_start'] or position >= align_attrs['mapped_end']:
        return None
    if len(raw_dat) > 99999999:
        raise ValueError(fast5_fn + ": max signal length exceed 99999999")
    if any(len(vals) <= 1 for vals in (corr_data, raw_dat)):
        raise NotImplementedError(('One or no segments or signal present in read.'))
    event_starts = corr_data['start'] + corr_start_rel_to_raw
    event_lengths = corr_data['length']
    event_bases = corr_data['base']
    fast5_data.close()
    # ~ label_data = np.array(
    # ~ list(zip(event_starts, event_lengths, event_bases)),
    # ~ dtype=[('start', '<u4'), ('length', '<u4'), ('base', 'S1')])
    return (raw_dat, event_bases, event_starts, event_lengths,start_position,strand,chrome)


def extract_feature(signal, event_start, event_length, base,start_position,end_position,strand,norm=True):
    global BASE_LIST
    position=FLAG.pos
    if strand == '+':
        current_index = position-start_position
    else:
        current_index = end_position - position - 1
    if current_index-FLAG.len < 0:
        start = 0
    else:
        start = current_index-FLAG.len
    if current_index + FLAG.len >= (end_position - start_position):
        end=len(event_start)
    else:
        end=current_index + FLAG.len + 1

    seq_array = np.array(list(base))
    assert seq_array.shape == event_length.shape
    # uniq_arr = np.unique(signal)
    # signal = (signal - np.median(uniq_arr)) / float(robust.mad(uniq_arr))
    # normalization

    signal = signal[event_start[0]:event_start[-1] + event_length[-1]]

    event_start = event_start-event_start[0]
    if norm:
        signal = normalize_signal_with_lim(signal)

    # filter too short or long dwell time
    # dwell_filter_pctls = (5, 95)
    # dwell_min, dwell_max = np.percentile(event_length, dwell_filter_pctls)
    # valid_bases = np.logical_and.reduce(
    #     (
    #         event_length > dwell_min,
    #         event_length < dwell_max,
    #         np.logical_not(np.isnan(event_length)),
    #     )
    # )
    # event_length[~valid_bases]=0


    total_feature_per_reads = []
    # if end >=event_start.shape[0]:
    #     #print(1)
    # draw_signal(norm_signal[event_start[start]:event_start[end]], event_start[start:end] - event_start[start],
    #             base[start:end])
    raw_signal_every = [signal[event_start[start + x]:event_start[start + x] + event_length[start + x]] for x in range(end-start)]

    for i, element in enumerate(raw_signal_every):
        if event_length[start + i] == 0:
            continue
        if strand == '+':
            final_position = position - FLAG.len + i
        else:
            final_position = position + FLAG.len - i
        temp = [np.mean(element), np.std(element), np.median(element), event_length[start + i], final_position]
        total_feature_per_reads.append(temp)
    total_feature_per_reads=np.array(total_feature_per_reads)
    return total_feature_per_reads

def extract_file(input_file,mode):
    basecall_group = FLAG.basecall_group
    basecall_subgroup = FLAG.basecall_subgroup

    chromosome = FLAG.chrom
    position = FLAG.pos
    strand = FLAG.strand
    try:
        temp_result = get_label_raw(input_file, basecall_group, basecall_subgroup,chromosome, position,strand)
        if temp_result is None:
            return None
        (raw_data, raw_label, raw_start, raw_length, start_position, strand,chrome) = temp_result
    except Exception as e:
        # print(str(e))
        return None

    # ~ print(input_file,raw_start,raw_length,raw_label)
    total_seq = "".join([x.decode() for x in raw_label])
    base_len=raw_label.shape[0]
    end_position = start_position+base_len
    if mode:
        raw_data = np.flip(raw_data)
    matrix_feature = extract_feature(raw_data, raw_start, raw_length, total_seq,start_position,end_position,strand,FLAG.norm)
    if matrix_feature is None:
        return None
    del raw_data

    return matrix_feature

def extract_group(args, total_fl,subsapmle_num=500):
    global FLAG,BASE_LIST
    BASE_LIST=None
    FLAG=args
    # feature_matrix = np.zeros((1000000,40))
    results_list = []
    ##########
    pool = multiprocessing.Pool(processes=int(args.cpu))
    ##########
    for fl in total_fl:
        result_per_read = pool.apply_async(extract_file, (fl,args.rna))
        results_list.append(result_per_read)
    pool.close()
    ############################
    pbar = tqdm(total=len(total_fl), position=0, leave=True)

    result_list=[]
    for result_per_read in results_list:
        temp= result_per_read.get()
        if temp is not None:
            feature_per_read = temp
            result_list.append(feature_per_read)

            del  feature_per_read
        pbar.update(1)
    #############################
    pool.join()
    pbar.close()

    num_aligned = len(result_list)
    if subsapmle_num < num_aligned:
        result_list=random.sample(result_list,subsapmle_num)
    final_feature=[]
    for item in result_list:
        final_feature.extend(item)
    df=pd.DataFrame(final_feature)


    if df.shape[0] == 0:
        raise Exception("can not find basecall_group or basecall_subgroup in fast5 files or there is no read aligned on the position")
    print('\nextracted ',num_aligned,' aligned reads from fast5 files')
    df.columns = ['Mean', 'STD', 'Median', 'Dwell time', 'position']

    # 转化为数值
    df['position'] = df['position'].astype(int).astype(str)

    if num_aligned > 50:
        dwell_filter_pctls = (0.5, 99.5)
        dwell_min, dwell_max = np.percentile(df['Dwell time'].values, dwell_filter_pctls)
        df = df[(df['Dwell time'] > dwell_min) & (df['Dwell time'] < dwell_max)]
        df.reset_index(inplace=True,drop=True)
    #     item_list = ['Mean', 'STD']
    #     for item in item_list:
    #         # collect data
    #         dwell_filter_pctls = (2.5, 97.5)
    #         dwell_min, dwell_max = np.percentile(df[item].values, dwell_filter_pctls)
    #         df = df[(df[item] > dwell_min) & (df[item] < dwell_max)]
    #         df.reset_index(inplace=True,drop=True)
    return df,num_aligned


def create_read_list_file(path,results_path):
    if path[-1] == '/':
        path = path[:-1]
    total_fl = []
    os.system("find " + path + " -name \"*.fast5\" >" + results_path + "/files.txt")
    read_name_list_file = results_path + '/files.txt'
    print("Generated fast5 list file: ", read_name_list_file)
    for i in open(read_name_list_file, "r"):
        total_fl.append(i.rstrip())
    cmd='rm '+read_name_list_file
    os.system(cmd)
    return total_fl


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     # parser.add_argument('--basecall_group', default="RawGenomeCorrected_000",
#     #                     help='The attribute group to extract the training data from. e.g. RawGenomeCorrected_000')
#     # parser.add_argument('--basecall_subgroup', default='BaseCalled_template',
#     #                     help='Basecall subgroup Nanoraw resquiggle into. Default is BaseCalled_template')
#     # parser.add_argument('-i',"--fast5", default='/data/Ecoli_23s/data/L_rep2/single',
#     #                     help="fast5_file")
#     # parser.add_argument('-c',"--control_fast5", default='/data/Ecoli_23s/data/IVT_negative/single',
#     #                     help="control_fast5_file")
#     # parser.add_argument('-o',"--output", default="/data/Ecoli_23s/tombo_results_2030_plus", help="output_file")
#     # parser.add_argument("--chrom", default='NR_103073.1',help="Gene or chromosome name(head of your fasta file)")
#     # parser.add_argument("--pos", default=2030, type=int,help="site of your interest")
#     # parser.add_argument("--len", default=10, type=int, help="region around the position")
#     # parser.add_argument("--strand", default="+", help="Strand of your interest")
#     # parser.add_argument('-t',"--cpu", default=4, type=int, help="num of process")
#     # parser.add_argument("--ref", default="/data/Ecoli_23s/23S_rRNA.fasta", help="fasta file")
#     # parser.add_argument("--overplot-number", default=500, type=int, help="Number of read will be used to plot")
#     parser.add_argument('--basecall_group', default="RawGenomeCorrected_000",
#                         help='The attribute group to extract the training data from. e.g. RawGenomeCorrected_000')
#     parser.add_argument('--basecall_subgroup', default='BaseCalled_template',
#                         help='Basecall subgroup Nanoraw resquiggle into. Default is BaseCalled_template')
#     parser.add_argument('-i',"--fast5", required=True,
#                         help="fast5_file")
#     parser.add_argument('-c',"--control_fast5",
#                         help="control_fast5_file")
#     parser.add_argument('-o',"--output", default="tombo_result", help="output_file")
#     parser.add_argument("--chrom", required=True, help="Gene or chromosome name(head of your fasta file)")
#     parser.add_argument("--pos", required=True, type=int,help="site of your interest")
#     parser.add_argument("--len", default=10, type=int, help="region around the position")
#     parser.add_argument("--strand", default="+", help="Strand of your interest")
#     parser.add_argument('-t',"--cpu", default=4, type=int, help="num of process")
#     parser.add_argument("--ref", required=True, help="fasta file")
#     parser.add_argument("--overplot-number", default=500, type=int, help="Number of read will be used to plot")
#     args = parser.parse_args()
#     args.pos = args.pos - 1
#     subsapmle_num = args.overplot_number
#     fasta=read_fasta_to_dic(args.ref)
#     base_list = fasta[args.chrom][args.pos-args.len:args.pos+args.len+1]
#     if args.strand == '-':
#         base_list="".join(list(reversed(base_list)))
#         base_list=reverse_fasta(base_list)
#     results_path = args.output
#     if not os.path.exists(results_path):
#         os.mkdir(results_path)
#
#     title = args.chrom + ':' + str(args.pos - args.len + 1) + '-' + str(args.pos + args.len + 2) + ':' + args.strand
#
#     wt_file=create_read_list_file(args.fast5,results_path)
#     df_wt,aligned_num_wt=extract_group(args, wt_file,subsapmle_num)
#     df_wt['type']='Sample'
#     try:
#         ivt_file=create_read_list_file(args.control_fast5,results_path)
#         df_ivt,aligned_num_ivt=extract_group(args, ivt_file,subsapmle_num)
#         df_ivt['type'] = 'Control'
#         df=pd.concat([df_wt,df_ivt])
#         title = title + '   Sample:' + str(aligned_num_wt) + '  Control:' + str(aligned_num_ivt)
#         category = pd.api.types.CategoricalDtype(categories=['Sample', "Control"], ordered=True)
#         df['type'] = df['type'].astype(category)
#     except:
#         args.control = None
#     if args.control_fast5 is None:
#         df = df_wt
#         df_wt['type'] = 'Single'
#         title = title + '   Sample:' + str(aligned_num_wt)
#
#     category_data = [str(args.pos + x) for x in range(-args.len, args.len + 1)]
#     category = pd.api.types.CategoricalDtype(categories=category_data, ordered=True)
#     df['position'] = df['position'].astype(category)
#
#     # draw_volin(df,results_path,args.pos,base_list,title)
#     # draw_boxplot(df,results_path,args.pos,base_list,title)
#     signal_plot(df, results_path, args.pos, base_list, title, 'merged')
#     signal_plot(df, results_path, args.pos, base_list, title,'boxplot')
#     signal_plot(df, results_path, args.pos, base_list, title, 'violin_plot')
#     print('\nsaved as ', args.output)
