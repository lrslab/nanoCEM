#!/usr/bin/env python
import pandas as pd
import os
import numpy as np
from tqdm import tqdm
import pysam
import argparse
import shutil
import time

def create_result_folder(output_path):
    print("Try to build generate the result path: "+output_path)
    if os.path.exists(output_path):
        print("Folder exists. It will be overwrite after 5 secs ...")
        time.sleep(5)
        print("Deleting folder...")
        shutil.rmtree(output_path)
        print("Folder deleted. Creating folder...")
        os.mkdir(output_path)
        print("Folder created.")
    else:
        print("Creating output folder...")
        os.mkdir(output_path)
        print("Folder created.")

def select_sub_reads(line):
    result = line[0].split("/")[-1].split(".")[0]
    return result

def copy2new_path(line,results_path):
    global pbar
    fold_name = line[0].split('/')[-2]
    if not os.path.exists(results_path+'/single/'+fold_name):
        os.mkdir(results_path+'/single/'+fold_name)
    os.system('cp ' + line[0] + ' ' + results_path+'/single/'+fold_name)
    pbar.update(1)

def transfer_fast5(total_fl, df,results_path):
    global pbar
    total_fl.columns=[0,'readname']
    df = df[['readname']]
    df[2] = '.'

    temp = pd.merge(total_fl, df, on='readname', how='inner')

    temp=temp.drop_duplicates(subset='readname')
    print("Start to extract "+str(temp.shape[0])+' reads ...')
    # temp=temp.sample(n=1544)
    pbar = tqdm(total=temp.shape[0], position=0, leave=True)


    os.mkdir(results_path+'/single/')
    temp.apply(copy2new_path,results_path=results_path, axis=1)
    pbar.close()

def main(args):
    fast5_path = args.fast5
    results_path = args.output
    sam_file=args.bam
    create_result_folder(results_path)
    print("Loading bam files ... ")
    if sam_file.split(".")[-1] == 'bam':
        sam_file=pysam.AlignmentFile(sam_file,'rb')
    else:
        sam_file=pysam.AlignmentFile(sam_file,'r')

    read_list=[]
    if args.chrom is None or args.pos is None:
        for read in sam_file.fetch():
            read_list.append(read.qname)
    else:
        for read in sam_file.fetch(args.chrom, args.pos - 10, args.pos + 10):
            read_list.append(read.qname)
    df = pd.DataFrame(read_list)
    df.columns=['readname']

    #READ FAST5 FILE LIST
    print("Collecting fast5 files ... ")
    os.system("find " + fast5_path + " -name \"*.fast5\" >" + results_path+"/files.txt")
    fast5_file =  results_path+"/files.txt"
    print("Generated fast5 list file: ", fast5_file)

    print("Subsampling fast5 files ... ")
    total_fl = []
    for i in open(fast5_file, "r"):
        total_fl.append(i.rstrip())
    total_fl = pd.DataFrame(np.array(total_fl))
    total_fl[1] = total_fl.apply(select_sub_reads, axis=1)
    transfer_fast5(total_fl, df, results_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--fast5", required=True,
                        help="fast5_path")
    parser.add_argument('-b',"--bam", required=True,
                        help="bam_path")
    parser.add_argument('-o',"--output", default='subsample_single', help="output_file")
    parser.add_argument("--chrom",help="Gene or chromosome name(head of your fasta file)")
    parser.add_argument("--pos", type=int, help="site of your interest")
    # parser.add_argument("-f","--fast5", default='/data/Ecoli_23s/data/L_rep2/single/',
    #                     help="fast5_path")
    # parser.add_argument('-b',"--bam", default="/data/Ecoli_23s/data/L_rep2/file.bam",
    #                     help="bam_path")
    # parser.add_argument('-o',"--output", default='/data/Ecoli_23s/data/L_rep2/subsample2', help="output_file")
    # parser.add_argument("--chrom", default='NR_103073.1',help="Gene or chromosome name(head of your fasta file)")
    # parser.add_argument("--pos", default=2029, type=int, help="site of your interest")
    args = parser.parse_args()
    main(args)
