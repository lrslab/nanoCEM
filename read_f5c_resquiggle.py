import re
import numpy as np
import pandas as pd
import pyslow5
import pysam
import plotnine as p9
from normalization import normalize_signal
import os
import argparse
from plot import draw_boxplot,draw_volin
def extract_feature(line,position):

    read_id = line[0]
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
    for item in moves:
        if 'D' in item:
            deletion = int(item[:-1])
            for i in range(deletion):
                event_length.append(0)
        elif 'I' in item:
            insertion = insertion + int(item[:-1])
            return None
        elif '=' in item:
            return None
        else:
            event_length.append(int(item))
    # build event_length from move table
    read = s5.get_read(read_id, aux=["read_number", "start_mux"], pA=True)
    start_index = line[2]
    end_index = line[3]
    event_length = np.array(event_length)
    #  assert len_raw_signal in paf and blow5
    try:
        assert end_index-start_index == np.sum(event_length)
        assert read['len_raw_signal'] == line[1]
    except Exception:
        print("Warning: 1 read's length of signal is not equal between blow5 and paf")
        return None

    # normalized signal
    signal = read['signal']
    signal=normalize_signal(signal)

    # create event start and flip all table
    signal = signal[start_index:end_index]
    signal = np.flip(signal)
    event_length = np.flip(event_length)
    event_start = event_length.cumsum()
    event_start = np.insert(event_start, 0, 0)[:-1]

    # index query and reference
    aligned_pair=info_dict[read_id]['pairs']
    qlen = info_dict[read_id]['query_length']
    gap = qlen - event_length.shape[0]
    aligned_pair[0] = aligned_pair[0] - gap
    aligned_pair = aligned_pair[aligned_pair[0] >= 0]
    if aligned_pair.shape[0]==0:
        return None
    read_pos = aligned_pair[0].values
    ref_pos = aligned_pair[1].values

    # extract raw signal by event length and event start
    total_feature_per_reads = []
    raw_signal_every = [signal[event_start[x]:event_start[x] + event_length[x]] for x in
                        read_pos]
    # calculate mean median and dwell time
    for i, element in enumerate(raw_signal_every):
        if event_length[read_pos[i]] == 0:
            continue
        if ref_pos[i] == position:
            if event_length[i] == 0:
                return None
        temp = [np.mean(element), np.std(element), np.median(element), event_length[read_pos[i]],ref_pos[i]]
        total_feature_per_reads.append(temp)
    return total_feature_per_reads

def extract_pairs_pos(bam_file,position,length,chromosome):

    result_dict={}
    for read in bam_file.fetch(chromosome,position-length,position+length+1):
        # if read.qname=='0eaa68b9-5989-44ca-8e00-52eba6ba3ccb':
        #     print(1)
        start_position=read.reference_start
        end_position=read.reference_end
        if position < start_position or position > end_position:
            continue
        # unit
        aligned_pair = np.array(read.aligned_pairs)
        aligned_pair = pd.DataFrame(aligned_pair)
        aligned_pair.dropna(inplace=True,ignore_index=True)
        aligned_pair = aligned_pair[(aligned_pair[1] >= position - length) & (aligned_pair[1] <= position + length)]


        temp={}
        temp['pairs']=aligned_pair
        temp['query_length'] = read.query_length
        result_dict[read.qname] = temp
    start = read.reference_start
    ref = read.get_reference_sequence()
    ref = ref[position-start-length:position-start+length+1]
    return result_dict,ref.upper()



def read_blow5(path,position,length,chromo=None,strand=None):
    global info_dict,s5
    bam_file=path+".bam"
    bam_file=pysam.AlignmentFile(bam_file,'rb')
    info_dict,ref=extract_pairs_pos(bam_file,position,length,chromo)

    slow5 = path+".blow5"
    s5 = pyslow5.Open(slow5, 'r')

    df=pd.read_csv(path+".paf",sep='\t',header=None)

    df["feature"] = df.apply(extract_feature,position=position,axis=1)
    df.dropna(inplace=True)
    num_aligned = df.shape[0]
    final_feature=[]
    for item in df["feature"]:
        final_feature.extend(item)
    final_feature=pd.DataFrame(final_feature)
    final_feature.columns=['Mean','STD','Median','Dwell_time','position']
    final_feature['position'] = final_feature['position'].astype(int).astype(str)
    return final_feature,num_aligned,ref




# plot = p9.ggplot(df, p9.aes(x='position', y=item)) \
#                + p9.geom_boxplot( outlier_shape='',size=0.25) \
#                + p9.theme_bw() \
#                + p9.theme(
#             figure_size=(6, 3),
#             panel_grid_minor=p9.element_blank(),
#             axis_text=p9.element_text(size=13),
#             axis_title=p9.element_text(size=13),
#             title=p9.element_text(size=13),
#             legend_position='none'
#         )
# print(plot)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", default='data/test_fast5/bhk/resquiggle_result/bhk',
                        help="blow5_path")
    parser.add_argument('-c',"--control", default='data/test_fast5/ivt/resquiggle_result/ivt',
                        help="control_blow5_path")
    parser.add_argument('-o',"--output", default="./f5c_results", help="output_file")
    parser.add_argument("--chrom", default='NC_001547.1',help="bed file to extract special site datasets")
    parser.add_argument("--pos", default=3929, help="bed file to extract special site datasets")
    parser.add_argument("--len", default=5, help="bed file to extract special site datasets")
    parser.add_argument("--strand", default="+", help="bed file to extract special site datasets")
    args = parser.parse_args()


    results_path = args.output
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    df_wt,aligned_num_wt,_=read_blow5(args.input,args.pos,args.len)

    df_wt['type']='Sample'
    df_ivt,aligned_num_ivt,base_list=read_blow5(args.control,args.pos,args.len)
    df_ivt['type']='Sample'
    df_ivt['type'] = 'Control'
    df=pd.concat([df_wt,df_ivt])

    category_data = [str(args.pos + x) for x in range(-args.len, args.len + 1)]
    category = pd.api.types.CategoricalDtype(categories=category_data, ordered=True)
    df['position'] = df['position'].astype(category)


    category = pd.api.types.CategoricalDtype(categories=['Sample',"Control"], ordered=True)
    df['type'] = df['type'].astype(category)
    df['Dwell_time'] = np.log10(df['Dwell_time'].values)

    draw_volin(df,results_path,args.pos,args.len,args.chrom,base_list,aligned_num_wt,aligned_num_ivt,strand="+")
    draw_boxplot(df,results_path,args.pos,args.len,args.chrom,base_list,aligned_num_wt,aligned_num_ivt,strand="+")
    print('\nsaved as ', args.output)
