import numpy as np
import h5py
import pandas as pd
import plotnine as p9
import argparse
import os
import multiprocessing
from tqdm import tqdm
from matplotlib import pyplot as plt
# from ._stats import c_new_mean_stds
from normalization import normalize_signal_with_lim
from plot import draw_volin,draw_boxplot
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial']


def matrix_append(total_matrix,new_matrix,index):
    len_new_matrix=new_matrix.shape[0]
    if index+len_new_matrix>total_matrix.shape[0]:
        temp_matrix=np.zeros((500000,5),dtype=float)
        total_matrix=np.vstack((total_matrix,temp_matrix))
    total_matrix[index:index+len_new_matrix]=new_matrix
    return total_matrix

def get_label_raw(fast5_fn, basecall_group, basecall_subgroup ,chromosome, position,strand):
    ##Open file
    # fast5_fn='/media/zhguo/QUO/ivt_positive/ivt_positive/single/42/878b9de7-043b-4d5e-b194-7771c5767ff5.fast5'
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


def extract_feature(signal, event_start, length, base,start_position,end_position):
    global BASE_LIST
    position=FLAG.pos
    current_index=position-start_position
    if current_index-FLAG.len < start_position:
        start = start_position
    else:
        start = current_index-FLAG.len
    if current_index + FLAG.len > end_position:
        end=end_position
    else:
        end=current_index + FLAG.len
    seq_array = np.array(list(base))
    assert seq_array.shape == length.shape
    # uniq_arr = np.unique(signal)
    # signal = (signal - np.median(uniq_arr)) / float(robust.mad(uniq_arr))
    # normalization
    norm_signal=normalize_signal_with_lim(signal)

    total_feature_per_reads = []
    base_list = base[start:end+1]

    raw_signal_every = [norm_signal[event_start[start + x]:event_start[start + x] + length[start + x]] for x in range(end-start+1)]

    for i, element in enumerate(raw_signal_every):
        if length[start + i] < 5:
            continue
        temp = [np.mean(element), np.std(element), np.median(element), length[start + i], position - FLAG.len +i]
        total_feature_per_reads.append(temp)
    total_feature_per_reads=np.array(total_feature_per_reads)
    return total_feature_per_reads,base_list

def extract_file(input_file):
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
    raw_data = raw_data[::-1]
    # ~ print(input_file,raw_start,raw_length,raw_label)
    total_seq = "".join([x.decode() for x in raw_label])
    base_len=raw_label.shape[0]
    end_position=start_position+base_len
    matrix_feature = extract_feature(raw_data, raw_start, raw_length, total_seq,start_position,end_position)
    if matrix_feature is None:
        return None
    del raw_data
    return matrix_feature

def extract_group(args, total_fl):
    global FLAG,BASE_LIST
    BASE_LIST=None
    FLAG=args
    output=args.output
    # feature_matrix = np.zeros((1000000,40))
    results_list = []
    ##########
    pool = multiprocessing.Pool(processes=int(args.cpu))
    ##########
    for fl in total_fl:
        result_per_read = pool.apply_async(extract_file, (fl,))
        results_list.append(result_per_read)
    pool.close()
    ############################
    pbar = tqdm(total=len(total_fl), position=0, leave=True)
    total_index=0
    feature_matrix = np.zeros((500000, 5),dtype=float)
    aligned_num=0
    for result_per_read in results_list:
        temp= result_per_read.get()
        if temp is not None:
            feature_per_read, base_list = temp
            aligned_num=aligned_num+1
            read_len=feature_per_read.shape[0]
            # npd_file.append(feature_per_read)
            feature_matrix=matrix_append(feature_matrix,feature_per_read,total_index)
            total_index=total_index+read_len

            del  feature_per_read
        pbar.update(1)
    #############################
    pool.join()
    pbar.close()
    if total_index == 0:
        raise Exception("can not find basecall_group or basecall_subgroup in fast5 files")
    print('\nextracted ',aligned_num,' aligned reads from fast5 files')

    feature_matrix= feature_matrix[:total_index]
    df=pd.DataFrame(feature_matrix)
    df.columns=['Mean','STD','Median','Dwell_time','position']
    # 去除非数字字符

    # 转化为数值
    df['position'] = df['position'].astype(int).astype(str)
    return df,base_list,aligned_num


def create_read_list_file(path,results_path):
    if path[0] == '/':
        path = path
    else:
        path = os.path.abspath('.') + '/' + path
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--basecall_group', default="RawGenomeCorrected_000",
                        help='The attribute group to extract the training data from. e.g. RawGenomeCorrected_000')
    parser.add_argument('--basecall_subgroup', default='BaseCalled_template',
                        help='Basecall subgroup Nanoraw resquiggle into. Default is BaseCalled_template')
    parser.add_argument('-i',"--fast5", default='data/test_fast5/bhk/single',
                        help="fast5_file")
    parser.add_argument('-c',"--control_fast5", default='data/test_fast5/ivt/single',
                        help="fast5_file")
    parser.add_argument('-o',"--output", default="./tombo_results", help="output_file")
    parser.add_argument("--chrom", default='NC_001547.1',help="bed file to extract special site datasets")
    parser.add_argument("--pos", default=3929, help="bed file to extract special site datasets")
    parser.add_argument("--len", default=5, help="bed file to extract special site datasets")
    parser.add_argument("--strand", default="+", help="bed file to extract special site datasets")
    parser.add_argument("--cpu", default=4, type=int, help="num of process")
    args = parser.parse_args()


    results_path = args.output
    if not os.path.exists(results_path):
        os.mkdir(results_path)
    wt_file=create_read_list_file(args.fast5,results_path)
    df_wt,_,aligned_num_wt=extract_group(args, wt_file)
    df_wt['type']='Sample'
    ivt_file=create_read_list_file(args.control_fast5,results_path)
    df_ivt,base_list,aligned_num_ivt=extract_group(args, ivt_file)
    df_ivt['type'] = 'Control'
    df=pd.concat([df_wt,df_ivt])

    category_data = [str(args.pos + x) for x in range(-args.len, args.len + 1)]
    category = pd.api.types.CategoricalDtype(categories=category_data, ordered=True)
    df['position'] = df['position'].astype(category)


    category = pd.api.types.CategoricalDtype(categories=['Sample',"Control"], ordered=True)
    df['type'] = df['type'].astype(category)
    df['Dwell_time'] = np.log10(df['Dwell_time'].values)


    draw_volin(df,results_path,args.pos,args.len,args.chrom,base_list,aligned_num_wt,aligned_num_ivt,strand="+")
    draw_volin(df,results_path,args.pos,args.len,args.chrom,base_list,aligned_num_wt,aligned_num_ivt,strand="+")
    print('\nsaved as ', args.output)
