from Bio import SeqIO
import os, sys
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        To get RSCU values from genbank files.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""

        return f"{description}\n\n{help_text}"


def RSCU_filter(input_file, filter_length=300):
    if filter_length is None:
        filter_length = 300
        print(f"The CDS filter length is set as default length: {filter_length} bp")
    else:
        print(f"The CDS filter length is {filter_length} bp")

    aa = []
    for rec in SeqIO.parse(input_file, 'genbank'):
        for feature in rec.features:
            if feature.type == 'CDS':
                if feature.qualifiers['gene'][0] not in aa:
                    aa.append(feature.qualifiers['gene'][0])
    print(f"Total gene number: {len(aa)}")

    # delete duplicated regions
    file_names, file_extension = os.path.splitext(os.path.basename(input_file))
    save_results_dir = file_names
    if os.path.exists(save_results_dir):
        print(f"The directory {save_results_dir} has been existed, please change the file name or delete the directory")
        sys.exit()
    os.mkdir(save_results_dir)
    remove_duplicated = os.path.join(save_results_dir, 'remove_duplicated.fasta')
    # print(remove_duplicated)
    with open(remove_duplicated, 'w') as ff:
        for i in aa:
            for rec in SeqIO.parse(input_file, format='genbank'):
                c = []
                for feature in rec.features:
                    if feature.type == 'CDS':
                        if feature.qualifiers['gene'][0] == i:
                            c.append(feature.extract(rec.seq))

                if len(c) == 1:
                    ff.write('>{}\n{}\n'.format(i, str(c[0])))
                if len(c) == 2:
                    c.remove(c[0]) if len(c[0]) <= len(c[1]) else c.remove(c[1])
                    ff.write('>{}\n{}\n'.format(i, str(c[0])))
    ff.close()
    # delete cds < 300 bp
    filter_length_file = os.path.join(save_results_dir, 'after_filter_' + str(filter_length) + '.fasta')
    filter_sequences = open(os.path.join(save_results_dir, 'filter_sequences.txt'), 'w')
    save_sequences = open(os.path.join(save_results_dir, 'save_sequences_name.txt'), 'w')
    with open(filter_length_file, 'w') as ff:
        number = 0
        for rec in SeqIO.parse(remove_duplicated, format='fasta'):
            if len(rec.seq) >= filter_length and rec.seq[:3] == 'ATG':
                save_sequences.write(f"{rec.id}\n")
                number += 1
                ff.write('>{}\n{}\n'.format(rec.id, str(rec.seq)))
            else:
                if len(rec.seq) < filter_length:
                    filter_sequences.write(f"{rec.id}\t seqeunce length samller than {filter_length} bp\n")
                if rec.seq[:3] != 'ATG':
                    filter_sequences.write(f"{rec.id}\t seqeunce not start with ATG\n")
        print(f"After fillter, only {number} sequences saved!")
        save_sequences.write(f"{number}")
        ff.close()
        save_sequences.close()
        filter_sequences.close()

    # merge sequences
    merge_fasta = os.path.join(save_results_dir, 'merge.fasta')
    merge_fasta_w = open(merge_fasta, 'w')
    merge_fasta_w.write(">all_merge_cds\n")
    for rec in SeqIO.parse(filter_length_file, format='fasta'):
        merge_fasta_w.write(str(rec.seq))
    merge_fasta_w.close()

    # final step: calculate rscu values
    codon_2_aa = {
        'TGG': 'Trp',   'GGT': 'Gly',   'AGG': 'Arg',   'CGA': 'Arg',
        'TGT': 'Cys',   'GGG': 'Gly',   'AGA': 'Arg',   'CGC': 'Arg',
        'TGC': 'Cys',   'GGA': 'Gly',   'AGC': 'Ser',   'CGG': 'Arg',
        'TAC': 'Tyr',   'GGC': 'Gly',   'AGT': 'Ser',   'CGT': 'Arg',
        'TAT': 'Tyr',   'GAT': 'Asp',   'AAG': 'Lys',   'CAA': 'Gin',
        'TCT': 'Ser',   'GAC': 'Asp',   'AAA': 'Lys',   'CAG': 'Gin',
        'TCG': 'Ser',   'GAA': 'Glu',   'AAC': 'Asn',   'CAT': 'His',
        'TCA': 'Ser',   'GAG': 'Glu',   'AAT': 'Asn',   'CAC': 'His',
        'TCC': 'Ser',   'GCA': 'Ala',   'ACC': 'Thr',   'CCT': 'Pro',
        'TTT': 'Phe',   'GCC': 'Ala',   'ACA': 'Thr',   'CCG': 'Pro',
        'TTC': 'Phe',   'GCG': 'Ala',   'ACG': 'Thr',   'CCA': 'Pro',
        'TTG': 'Leu',   'GCT': 'Ala',   'ACT': 'Thr',   'CCC': 'Pro',
        'TTA': 'Leu',   'GTA': 'Val',   'ATC': 'Ile',   'CTA': 'Leu',
        'TAG': 'Ter',   'GTC': 'Val',   'ATA': 'Ile',   'CTT': 'Leu',
        'TGA': 'Ter',   'GTG': 'Val',   'ATT': 'Ile',   'CTC': 'Leu',
        'TAA': 'Ter',   'GTT': 'Val',   'ATG': 'Met',   'CTG': 'Leu',

        }
    seqs = [rec.seq for rec in SeqIO.parse(merge_fasta, format='fasta')][0]
    print(f"After merging, the total CDS length is  {len(seqs)}")
    aa_list = []
    seqs = str(seqs)
    # stat number of codons
    for i in range(0, len(seqs), 3):
        amino_codon = seqs[i:i + 3]
        aa_list.append(amino_codon)
    number_codon = []
    for key, value in codon_2_aa.items():
        # print(key, aa_list.count(key), value)
        number_codon.append([key, aa_list.count(key), value])
    # print(number_codon)
    # calculate RSCU values
    save_results_file = os.path.join(os.path.dirname(merge_fasta), 'RSCU_results.txt')
    calculate_results = open(save_results_file, 'w')
    all_single_amino_number = {}
    all_same_number = {}
    for i in range(len(number_codon)):
        amino_number = number_codon[i][1]
        amino_name = number_codon[i][2]
        if amino_name not in all_single_amino_number.keys():
            all_single_amino_number[amino_name] = amino_number
            all_same_number[amino_name] = 1
        else:
            all_single_amino_number[amino_name] += amino_number
            all_same_number[amino_name] += 1
    for i in range(len(number_codon)):
        codon_name = number_codon[i][0]
        codon_number = number_codon[i][1]
        amino_name = number_codon[i][2]
        # rfsc_value = round(codon_number / (all_single_amino_number[amino_name]), 2)
        rscu_value = round(codon_number/(all_single_amino_number[amino_name]/all_same_number[amino_name]), 2)
        calculate_results.write(f"{amino_name}\t{codon_name}\t{codon_number}\t{rscu_value}\n")
        # print(amino_name, codon_name, codon_number, rfsc_value)
    calculate_results.close()

    ## the introduction for each directory
    info_file =os.path.join(save_results_dir, 'introduction.txt')
    with open(info_file, 'w') as ff:
        ff.write(
            'Introduction to each file:\n'
            'remove_duplicated.fasta:\t Deleted duplicated genes and their sequences.\n'
            'filter_sequences.txt:\t Deleted sequences starting with non-ATG and shorter than the specified length.\n'
            'after_filter_300.fasta:\t All remaining sequences after filtering.\n'
            'save_sequences_name.txt:\t Names and numbers of all remaining genes after filtering.\n'
            'merge.fasta:\t Merged sequences of all remaining gene sequences.\n'
            'RSCU_results.txt:\t Number of each codon and RSCU values. \n'
            '\tEach row has four columns: amino acid, corresponding codon, number of codons, and RSCU value.'
        )





if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-d", "--work_dir", help="Input directory of genbank files")
    parser.add_argument('-l', '--filter_length', help='CDS filter length, default is 300', type=int)
    args = parser.parse_args()
    if args.work_dir:
        for i in os.listdir(args.work_dir):
            if i.endswith('gb') or i.endswith('gbk'):
                file_path = os.path.abspath(os.path.join(args.work_dir, i))
                print('-'*80)
                print(file_path)
                RSCU_filter(file_path, args.filter_length)
                print('-'*80)
            else:
                print(f"please input genbank format files and endwith 'gb' or 'gbk'! ")
    else:
        parser.print_help()
