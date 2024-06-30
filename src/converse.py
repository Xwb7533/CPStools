from Bio import SeqIO
import os
import sys
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Converse genbank format files to fasta/tbl/mVISTA format.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""
        return f"{description}\n\n{help_text}"


def parse_genbank_file(input_file):
    try:
        with open(input_file, 'r') as f:
            records = list(SeqIO.parse(f, 'genbank'))
            if records:
                return records
            else:
                print(f"No valid GenBank records found in {input_file}")
                return None
    except (ValueError, IndexError, KeyError, IOError, FileNotFoundError) as e:
        print(f"Failed to parse GenBank file '{input_file}': {e}")
        return None


def gb2fa(input_file, save_file):
    if parse_genbank_file(input_file):
        filename_with_extension = os.path.basename(input_file)
        seq_id, _ = os.path.splitext(filename_with_extension)
        records = parse_genbank_file(input_file)
        with open(save_file, 'w') as ff:
            for rec in SeqIO.parse(input_file, format='genbank'):
                ff.write(f'>{seq_id}\n{rec.seq}\n')
            ff.close()
    else:
        sys.exit()


def gb2mVISTA(input_file, save_file):


    all_info = []
    if parse_genbank_file(input_file):
        for rec in SeqIO.parse(input_file, format='genbank'):
        # records = parse_genbank_file(input_file)
            for feature in rec.features:
                if feature.type == 'gene':
                    for part in feature.location.parts:
                        if int(part.strand) == 1:
                            if ['>', int(part.start) + 1, int(part.end), feature.qualifiers['gene'][0]] not in all_info:
                                all_info.append(['>', int(part.start) + 1, int(part.end), feature.qualifiers['gene'][0]])
                        else:
                            if ['<', int(part.start) + 1, int(part.end), feature.qualifiers['gene'][0]] not in all_info:
                                all_info.append(['<', int(part.start) + 1, int(part.end), feature.qualifiers['gene'][0]])
                elif feature.type == 'CDS':
                    for part in feature.location.parts:
                        if int(part.strand) == 1:
                            if [int(part.start) + 1, int(part.end), 'exon'] not in all_info:
                                all_info.append([int(part.start) + 1, int(part.end), 'exon'])
                        else:
                            if [int(part.start) + 1, int(part.end), 'exon'] not in all_info:
                                all_info.append([int(part.start) + 1, int(part.end), 'exon'])
                elif feature.type == 'rRNA' or feature.type == 'tRNA':
                    for part in feature.location.parts:
                        if int(part.strand) == 1:
                            all_info.append([int(part.start) + 1, int(part.end), 'utr'])
                        else:
                            all_info.append([int(part.start) + 1, int(part.end), 'utr'])
                else:
                    pass
    with open(save_file, 'w') as ff:
        for list_ in all_info:
            if len(list_) == 4:
                ff.write(f"{list_[0]} {list_[1]} {list_[2]} {list_[3]}\n")
            else:
                ff.write(f"{list_[0]} {list_[1]} {list_[2]}\n")


def gb2tbl(input_file, save_file):
    output_order = ["trans_splicing", "exception", "pesudo", "codon_start", "product", "gene", "transl_table"]
    if parse_genbank_file(input_file):
        with open(save_file, 'w') as ff:
            for rec in SeqIO.parse(input_file, 'genbank'):
                ff.write(f">Feature {rec.id}\n")
                for feature in rec.features[1:]:
                    if len(feature.location.parts) == 1:
                        for part in feature.location.parts:
                            if part.strand == 1:
                                ff.write(f"{part.start + 1}\t{part.end}\t{feature.type}\n")
                                for key in output_order:
                                    if key not in feature.qualifiers:
                                        continue
                                    value = feature.qualifiers[key][0]
                                    ff.write(f"\t\t\t{key}\t{value}\n")
                            else:
                                ff.write(f"{part.end}\t{part.start + 1}\t{feature.type}\n")
                                for key in output_order:
                                    if key not in feature.qualifiers:
                                        continue
                                    value = feature.qualifiers[key][0]
                                    ff.write(f"\t\t\t{key}\t{value}\n")
                    elif len(feature.location.parts) == 2:
                        if feature.location.strand == -1:
                            ff.write(f"{feature.location.parts[0].end}\t{feature.location.parts[0].start + 1}\t{feature.type}\n"
                                  f"{feature.location.parts[1].end}\t{feature.location.parts[1].start + 1}\n")
                            for key in output_order:
                                if key not in feature.qualifiers:
                                    continue
                                value = feature.qualifiers[key][0]
                                ff.write(f"\t\t\t{key}\t{value}\n")
                        elif feature.location.strand == 1:
                            ff.write(f"{feature.location.parts[0].start + 1}\t{feature.location.parts[0].end}\t{feature.type}\n"
                                  f"{feature.location.parts[1].start + 1}\t{feature.location.parts[1].end}\n")
                            for key in output_order:
                                if key not in feature.qualifiers:
                                    continue
                                value = feature.qualifiers[key][0]
                                ff.write(f"\t\t\t{key}\t{value}\n")
                        else:
                            ff.write(f"{feature.location.parts[0].start + 1}\t{feature.location.parts[0].end}\t{feature.type}\n"
                                  if feature.location.parts[0].strand == 1 else
                                  f"{feature.location.parts[0].end}\t{feature.location.parts[0].start + 1}\t{feature.type}\n")
                            ff.write(f"{feature.location.parts[1].start + 1}\t{feature.location.parts[1].end}\n"
                                  if feature.location.parts[1].strand == 1 else
                                  f"{feature.location.parts[1].end}\t{feature.location.parts[1].start + 1}\n")
                            for key in output_order:
                                if key not in feature.qualifiers:
                                    continue
                                value = feature.qualifiers[key][0]
                                ff.write(f"\t\t\t{key}\t{value}\n")
                    else:
                        if feature.location.strand == -1:
                            ff.write(f"{feature.location.parts[0].end}\t{feature.location.parts[0].start + 1}\t{feature.type}\n"
                                  f"{feature.location.parts[1].end}\t{feature.location.parts[1].start + 1}\n"
                                  f"{feature.location.parts[2].end}\t{feature.location.parts[2].start + 1}\n")
                            for key in output_order:
                                if key not in feature.qualifiers:
                                    continue
                                value = feature.qualifiers[key][0]
                                ff.write(f"\t\t\t{key}\t{value}\n")
                        elif feature.location.strand == 1:
                            ff.write(f"{feature.location.parts[0].start + 1}\t{feature.location.parts[0].end}\t{feature.type}\n"
                                  f"{feature.location.parts[1].start + 1}\t{feature.location.parts[1].end}\n"
                                  f"{feature.location.parts[2].start + 1}\t{feature.location.parts[2].end}\n")
                            for key in output_order:
                                if key not in feature.qualifiers:
                                    continue
                                value = feature.qualifiers[key][0]
                                ff.write(f"\t\t\t{key}\t{value}\n")
                        else:
                            ff.write(f"{feature.location.parts[0].start + 1}\t{feature.location.parts[0].end}\t{feature.type}\n"
                                  if feature.location.parts[0].strand == 1 else
                                  f"{feature.location.parts[0].end}\t{feature.location.parts[0].start + 1}\t{feature.type}\n")
                            ff.write(f"{feature.location.parts[1].start + 1}\t{feature.location.parts[1].end}\n"
                                  if feature.location.parts[1].strand == 1 else
                                  f"{feature.location.parts[1].end}\t{feature.location.parts[1].start + 1}\n")
                            ff.write(f"{feature.location.parts[2].start + 1}\t{feature.location.parts[2].end}\n"
                                  if feature.location.parts[2].strand == 1 else
                                  f"{feature.location.parts[2].end}\t{feature.location.parts[2].start + 1}\n")
                            for key in output_order:
                                if key not in feature.qualifiers:
                                    continue
                                value = feature.qualifiers[key][0]
                                ff.write(f"\t\t\t{key}\t{value}\n")
            ff.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-d", "--input_dir", help="Input path of genbank file", required=True)
    parser.add_argument("-m", "--mode", choices=['fasta', 'mVISTA', 'tbl'],
                        help="Mode: fasta for converting genbank format file into fasta format file;\n"
                             "mVISTA for converting genbank format file into mVISTA format file;\n"
                             "tbl for converting genbank format file into tbl format file", required=True)
    args = parser.parse_args()

    if args.input_dir:
        mode = args.mode
        output_dir = os.path.join(os.path.dirname(args.input_dir), mode)
        abs_output_dir = os.path.abspath(output_dir)
        
        if not os.path.exists(abs_output_dir):
            os.makedirs(abs_output_dir)
        else:
            print(f"{abs_output_dir} has existed, please select another one.")
            sys.exit()
        
        for file_ in os.listdir(args.input_dir):
            file_path = os.path.join(args.input_dir, file_)
            if file_path.endswith('gb') or file_path.endswith('gbk'):
                file_name, file_extension = os.path.splitext(os.path.basename(file_path))
                if mode == 'fasta':
                    save_file = os.path.join(abs_output_dir, file_name + '.fasta')
                    gb2fa(file_path, save_file)
                    print(f"{os.path.abspath(file_path)} has been converted into fasta format!")
                elif mode == 'mVISTA':
                    save_file = os.path.join(abs_output_dir, file_name + '.mVISTA')
                    gb2mVISTA(file_path, save_file)
                    print(f"{os.path.abspath(file_path)} has been converted into mVISTA format!")
                elif mode == 'tbl':
                    save_file = os.path.join(abs_output_dir, file_name + '.tbl')
                    gb2tbl(file_path, save_file)
                    print(f"{os.path.abspath(file_path)} has been converted into tbl format!")
                else:
                    print("Three modes are provided, you must specify one of them")
            else:
                print(f"{file_path} is not in genbank format, Please ensure it ends with 'gb' or 'gbk'")
    else:
        parser.print_help()
    