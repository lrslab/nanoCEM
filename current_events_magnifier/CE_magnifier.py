#!/usr/bin/env python
import argparse
import os
from current_events_magnifier.cem_utils import read_fasta_to_dic,reverse_fasta
import pandas as pd
from current_events_magnifier.plot import signal_plot
def init_parser():
    def add_public_argument(parser_input):
        parser_input.add_argument("--chrom", required=True, help="Gene or chromosome name(head of your fasta file)")
        parser_input.add_argument("--pos", required=True, type=int, help="site of your interest")
        parser_input.add_argument("--len", default=10, type=int, help="region around the position")
        parser_input.add_argument("--strand", default="+", help="Strand of your interest")
        parser_input.add_argument("--ref", required=True, help="fasta file")
        parser_input.add_argument("--overplot-number", default=500, type=int,
                                  help="Number of read will be used to plot")
    # Define the argument parser
    parser = argparse.ArgumentParser(description='A sample tool designed to visualize the features that distinguish between two groups of ONT data at the site level. It supports two re-squiggle pipeline(Tombo and f5c).')
    subparsers = parser.add_subparsers(dest='function')

    # tombo subparser
    parser_tombo = subparsers.add_parser('tombo', help='tackle tombo re-squiggle')
    parser_tombo.add_argument('--basecall_group', default="RawGenomeCorrected_000",
                        help='The attribute group to extract the training data from. e.g. RawGenomeCorrected_000')
    parser_tombo.add_argument('--basecall_subgroup', default='BaseCalled_template',
                        help='Basecall subgroup Nanoraw resquiggle into. Default is BaseCalled_template')
    parser_tombo.add_argument('-i', "--input_fast5", required=True,
                        help="input_fast5_file")
    parser_tombo.add_argument('-c', "--control_fast5",
                        help="control_fast5_file")
    parser_tombo.add_argument('-o', "--output", default="tombo_result", help="output_file")
    parser_tombo.add_argument('-t', "--cpu", default=4, type=int, help="num of process")
    add_public_argument(parser_tombo)

    # f5c subparser
    parser_f5c = subparsers.add_parser('f5c', help='tackle f5c re-squiggle')
    parser_f5c.add_argument("-i", "--input", required=True,
                        help="blow5_path")
    parser_f5c.add_argument('-c', "--control", required=True,
                        help="control_blow5_path")
    parser_f5c.add_argument('-o', "--output", default="f5c_result", help="output_file")
    add_public_argument(parser_f5c)
    return parser

if __name__ == '__main__':
    # Parse the arguments
    parser = init_parser()
    args = parser.parse_args()

    args.pos = args.pos - 1
    subsample_num = args.overplot_number
    fasta = read_fasta_to_dic(args.ref)
    base_list = fasta[args.chrom][args.pos - args.len:args.pos + args.len + 1]
    if args.strand == '-':
        base_list = "".join(list(reversed(base_list)))
        base_list = reverse_fasta(base_list)
    results_path = args.output
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    title = args.chrom + ':' + str(args.pos - args.len + 1) + '-' + str(args.pos + args.len + 2) + ':' + args.strand
    if args.function == 'tombo' :

        from  current_events_magnifier.read_tombo_resquiggle import create_read_list_file,extract_group
        wt_file = create_read_list_file(args.input_fast5, results_path)
        df_wt, aligned_num_wt = extract_group(args, wt_file, subsample_num)
        df_wt['type'] = 'Sample'
        try:
            ivt_file = create_read_list_file(args.control_fast5, results_path)
            df_ivt, aligned_num_ivt = extract_group(args, ivt_file, subsample_num)
            df_ivt['type'] = 'Control'
            df = pd.concat([df_wt, df_ivt])
            title = title + '   Sample:' + str(aligned_num_wt) + '  Control:' + str(aligned_num_ivt)
            category = pd.api.types.CategoricalDtype(categories=['Sample', "Control"], ordered=True)
            df['type'] = df['type'].astype(category)
        except:
            args.control_fast5 = None
        if args.control_fast5 is None:
            df = df_wt
            df_wt['type'] = 'Single'
            title = title + '   Sample:' + str(aligned_num_wt)
    elif args.function == 'f5c':
        from current_events_magnifier.read_f5c_resquiggle import read_blow5
        df_wt, aligned_num_wt = read_blow5(args.input, args.pos, args.len, args.chrom, args.strand, subsample_num)
        df_wt['type'] = 'Sample'
        try:
            df_ivt, aligned_num_ivt = read_blow5(args.control, args.pos, args.len, args.chrom, args.strand,
                                                 subsample_num)
            df_ivt['type'] = 'Control'

            df = pd.concat([df_wt, df_ivt])
            category = pd.api.types.CategoricalDtype(categories=['Sample', "Control"], ordered=True)
            df['type'] = df['type'].astype(category)

            title = title + '   Sample:' + str(aligned_num_wt) + '  Control:' + str(aligned_num_ivt)
        except:
            args.control = None
        if args.control is None:
            df = df_wt
            df_wt['type'] = 'Single'
            title = title + '   Sample:' + str(aligned_num_wt)

    category_data = [str(args.pos + x) for x in range(-args.len, args.len + 1)]
    category = pd.api.types.CategoricalDtype(categories=category_data, ordered=True)
    df['position'] = df['position'].astype(category)

    # draw_volin(df,results_path,args.pos,base_list,title)
    # draw_boxplot(df,results_path,args.pos,base_list,title)
    signal_plot(df, results_path, args.pos, base_list, title, 'merged')
    signal_plot(df, results_path, args.pos, base_list, title, 'boxplot')
    signal_plot(df, results_path, args.pos, base_list, title, 'violin_plot')
    print('\nsaved as ', args.output)
