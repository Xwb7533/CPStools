import os
import re
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = '''Adjust the LSC and SSC start in chloroplast genomes.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com'''

        return f"{description}\n\n{help_text}"

def adjust_SSC_forward(work_dir, save_dir, info_file):
    info_cont = open(info_file, 'r')
    info_line = info_cont.readline()
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    else:
        print(f"Directory '{save_dir}' already exists, please change another output directory!")
        sys.exit()
    while info_line:
        file_name = info_line.split('\t')[0]
        SSC_loc_start = int(info_line.split('\t')[3].split('-')[0].split(':')[1]) - 1
        SSC_loc_end = int(info_line.split('\t')[3].split('-')[1].strip())
        print(file_name)
        work_file = os.path.join(work_dir, file_name)
        save_file = os.path.join(save_dir, file_name)
        if os.path.exists(work_file):
            with open(save_file, 'w') as ff:
                for rec in SeqIO.parse(work_file, 'fasta'):
                    SSC_seq = Seq(str(rec.seq[SSC_loc_start:SSC_loc_end])).reverse_complement()
                    rev_seq = str(rec.seq[:SSC_loc_start]) + str(SSC_seq) + str(rec.seq[SSC_loc_end:])
                    ff.write(f'>{rec.id}\n{rev_seq}\n')
                info_line = info_cont.readline()
        else:
            abs_path = os.path.abspath(work_file)
            print(f"No such file: {abs_path}")
            info_line = info_cont.readline()
    print("Done")


def adjust_start_to_LSC(work_dir, save_dir, info_text):
    info_cont = open(info_text, 'r')
    info_line = info_cont.readline()
    # check save directory exist or not
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    else:
        print(f"Directory '{save_dir}' already exists, please change another output directory!")
        sys.exit()
    while info_line:
        file_name = info_line.split('\t')[0]
        LSC_loc = info_line.split('\t')[1].split('-')[0].split(':')[1]

        print(file_name)
        save_name = file_name

        work_file = os.path.join(work_dir, file_name)
        save_file = os.path.join(save_dir, save_name)

        if os.path.exists(work_file):
            with open(save_file, 'w') as ff:
                for rec in SeqIO.parse(work_file, 'fasta'):
                    if LSC_loc == '1':
                        ff.write(f'>{file_name}\n{rec.seq}')
                        info_line = info_cont.readline()
                    else:
                        adj_seq = str(rec.seq[int(LSC_loc) - 1:]) + str(rec.seq[:int(LSC_loc) - 1])
                        ff.write(f'>{file_name}\n{adj_seq}\n')
                        info_line = info_cont.readline()
        else:
            abs_path = os.path.abspath(work_file)
            print(f"No such file: {abs_path}")
            info_line = info_cont.readline()
    print("Done")

def adjust_start_to_RP(work_dir, save_dir, info_text):
    info_cont = open(info_text, 'r')
    info_line = info_cont.readline()
        # check save directory exist or not
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    else:
        print(f"Directory '{save_dir}' already exists, please change another output directory!")
        sys.exit()
    while info_line:

        def find_max_region_number(info_line):
            region_number_pairs = re.findall(r':(\d+)-(\d+)', info_line)
            max_region_number = max(max(int(num) for num in pair) for pair in region_number_pairs)
            return max_region_number


        seq_length = find_max_region_number(info_line)
        file_name = info_line.split('\t')[0]
        if len(info_line.split('\t')[1].split('-')) == 2:
            LSC_end = info_line.split('\t')[1].split('-')[1]
        else:
            LSC_end = info_line.split('\t')[1].split('-')[2]

        print(file_name)
        save_name = file_name

        work_file = os.path.join(work_dir, file_name)
        save_file = os.path.join(save_dir, save_name)

        if os.path.exists(work_file):
            with open(save_file, 'w') as ff:
                for rec in SeqIO.parse(work_file, 'fasta'):
                    my_seq = rec.seq.reverse_complement()
                    rp_seq = my_seq[seq_length-int(LSC_end):] + my_seq[:seq_length-int(LSC_end)]
                    ff.write(f'>{file_name}\n{rp_seq}\n')
                    info_line = info_cont.readline()
        else:
            abs_path = os.path.abspath(work_file)
            print(f"No such file: {abs_path}")
            info_line = info_cont.readline()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--work_dir", help="Input directory of fasta files")
    parser.add_argument("-o", "--save_dir", help="Output directory for save files")
    parser.add_argument("-f", "--info_file", help="file of the information of the fasta file four regions")
    parser.add_argument("-m", "--mode", choices=['SSC', 'LSC', 'RP'],
        help="Mode: SSC for adjust_SSC_forward, LSC for adjust_start_to_LSC, RP for adjust sequence to reverse_complement")
    args = parser.parse_args()

    if args.work_dir and args.save_dir and args.info_file:
        # Assuming the user specifies which function to call via an additional command line argument
        if args.mode == 'SSC':
            adjust_SSC_forward(args.work_dir, args.save_dir, args.info_file)
        elif args.mode == 'LSC':
            adjust_start_to_LSC(args.work_dir, args.save_dir, args.info_file)
        elif args.mode == 'RP':
            adjust_start_to_RP(args.work_dir, args.save_dir, args.info_file)
        else:
            print("Three modes are provided, you must specify one of them")
    else:
        parser.print_help()
