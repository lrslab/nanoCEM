from Bio import SeqIO
import pysam
import argparse

def main(args):
    result_dict= {}
    # read bam
    sam_file=args.bam
    print("Loading bam files ... ")
    if sam_file.split(".")[-1] == 'bam':
        sam_file = pysam.AlignmentFile(sam_file, 'rb')
    else:
        sam_file = pysam.AlignmentFile(sam_file, 'r')

    if args.chrom is None or args.pos is None:
        for read in sam_file.fetch():
            result_dict[read.qname]=1
    else:
        for read in sam_file.fetch(args.chrom, args.pos - 10, args.pos + 10):
            result_dict[read.qname]=1


    fastq_path =args.fastq

    # 创建一个新的FASTQ文件并写入序列记录
    if args.output.split(".")[-1] != 'fastq' and args.output.split(".")[-1]!='fq':
        print(args.output.split(".")[-1])
        args.output=args.output+'.fastq'
    new_fastq= open(args.output, "w")
    result_record_list=[]
    print("Loading fastq file ... ")
    with open(fastq_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # 处理每个序列记录
            if record.id in result_dict:
                result_record_list.append(record)
    print("Generating fastq file ... ")
    SeqIO.write(result_record_list, new_fastq, "fastq")
    new_fastq.close()
    print("Finished")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--fastq", required=True,
                        help="fast5_path")
    parser.add_argument('-b',"--bam", required=True,
                        help="bam_path")
    parser.add_argument('-o',"--output", default='subsample_single', help="output_file")
    parser.add_argument("--chrom",help="Gene or chromosome name(head of your fasta file)")
    parser.add_argument("--pos", type=int, help="site of your interest")
    args = parser.parse_args()
    main(args)