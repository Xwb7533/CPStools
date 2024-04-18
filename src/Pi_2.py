import os
import sys 
import argparse
from Bio import SeqIO


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Sorted intergenic spaces results as cp genome order.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""
        return f"{description}\n\n{help_text}"


def calculate_Pi_values(work_dir):
    if not work_dir.endswith('/'):
        work_dir += '/'
    all_pi_results = []
    for align_fasta in os.listdir(work_dir):
        if align_fasta.endswith('.fasta'):
            a, b, c, d = [], [], [], []
            pi = 0
            fasta_file = os.path.join(work_dir, align_fasta)
            for rec in SeqIO.parse(fasta_file, format='fasta'):
                a.append(rec.id)
                for i in range(len(rec.seq)):
                    if rec.seq[i] == '-':
                        b.append(i)
            all_number = len(a) * (len(a) - 1) / 2
            # delete all have '-' in seq location
            all_del = sorted(set(b))[::-1]
            for rec in SeqIO.parse(fasta_file, 'fasta'):
                for x in all_del:
                    rec.seq = list(rec.seq)
                    del rec.seq[x]
                d.append(rec.seq)
            # statistics same and diff
            for y in range(len(d[0])):
                c = []
                for x in d:
                    c.append(x[y])
                diff = 0
                for sig in range(len(a)):
                    for sig2 in range(sig + 1, len(a)):
                        if c[sig] != c[sig2]:
                            diff += 1
                pi += diff / all_number
            if len(d[0]) == 0:
                final_pi = 0
            else:
                final_pi = format(pi / len(d[0]), '.5f')
            all_pi_results.append(f"{align_fasta[:-6]}\t{final_pi}")
            print(f"{align_fasta[:-6]}\t{final_pi}")
    pi_results = os.path.join(work_dir, 'Pi_results.txt')
    with open(pi_results, 'w') as ff:
        for each_pi in all_pi_results:
            ff.write(f"{each_pi}\n")
    return pi_results


def IGS_sort_as_cp_order(input_file1, input_file2):
    pi_results = open(input_file1, 'r')
    cp_order_results = open(input_file2, 'r')
    results_file_path = os.path.join(os.path.dirname(input_file1), 'IGS_sort_as_cp_order.txt')
    reuslts_file = open(results_file_path, 'w')
    file1_line_list = pi_results.readlines()
    file2_line_list = cp_order_results.readlines()
    for IGS2 in file2_line_list:
        for IGS1 in file1_line_list:
            if IGS1.split('\t')[0] == IGS2.strip():
                print(IGS1, end='')
                reuslts_file.write(IGS1)
    reuslts_file.close()

def gene_sort_as_cp_order(input_file1, input_file2):
    pi_results = open(input_file1, 'r')
    cp_order_results = open(input_file2, 'r')
    results_file_path = os.path.join(os.path.dirname(input_file1), 'gene_sort_as_cp_order.txt')
    reuslts_file = open(results_file_path, 'w')
    file1_line_list = pi_results.readlines()
    file2_line_list = cp_order_results.readlines()
    for IGS2 in file2_line_list:
        for IGS1 in file1_line_list:
            if IGS1.split('\t')[0] == IGS2.strip():
                print(IGS1, end='')
                reuslts_file.write(IGS1)
    reuslts_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--input', type=str, help="Input the directory path of the multi-alignment sequences")
    parser.add_argument('-r', '--reference', type=str, help="Input the file path of 'cp_sort_IGS.txt'")
    parser.add_argument('-m', '--mode', type=str, choices=['IGS', 'gene'], 
        help="Mode: IGS for sorting intergenic spaces sequences into cp order;\n"
        "gene for sorting gene sequences into cp order;\n")
    args = parser.parse_args()
    print(args)
    if not (args.input and args.reference and args.mode):
        parser.print_help()
        sys.exit()
    else:
        pi_results_path = calculate_Pi_values(args.input)
        if args.mode == 'IGS':
            IGS_sort_as_cp_order(pi_results_path, args.reference)
            work_dir = os.path.dirname(pi_results_path)

            final_path = os.path.abspath(os.path.join(work_dir, "IGS_sort_as_cp_order.txt"))
            print(f"The sorted results have been written into:\n"
                f"\t\t\t{final_path}")
        elif args.mode == 'gene':
            gene_sort_as_cp_order(pi_results_path, args.reference)
            work_dir = os.path.dirname(pi_results_path)
            
            final_path = os.path.abspath(os.path.join(work_dir, "gene_sort_as_cp_order.txt"))
            print(f"The sorted results have been written into:\n"
                f"\t\t\t{final_path}")
        else:
            parser.print_help()
            sys.exit()
        


