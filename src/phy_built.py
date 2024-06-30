from Bio import SeqIO
import os
import sys
import argparse

class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Extract and sort common cds/protein sequences for phylogenetic analysis.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""
        return f"{description}\n\n{help_text}"


def find_ref(input_dir):
    # 找到目录中第一个以.gb或.gbk结尾的文件作为参考文件
    for filename in os.listdir(input_dir):
        if filename.endswith('.gb') or filename.endswith('.gbk'):
            ref_file = os.path.join(input_dir, filename)
            return ref_file
    
    # 如果没有找到符合条件的文件，返回 None 并打印提示信息
    print("No GenBank file found in the directory.")
    return None


def pro_for_phy(ref_file, output_dir):
    a = []
    work_dir = os.path.dirname(ref_file)
    for rec in SeqIO.parse(ref_file, format='genbank'):
        for feature in rec.features:
            if feature.type == 'CDS':
                if feature.qualifiers.get('gene') and feature.qualifiers.get('translation'):
                    if feature.qualifiers['gene'][0].lower() not in a:
                        a.append(feature.qualifiers['gene'][0].lower())
                else:
                    print(f"{feature} in genbank file has no gene name or translation!, please verify!!!")
    # delete not common cds
    for i in os.listdir(work_dir):
        file_path = os.path.join(work_dir, i)
        b = []
        for rec in SeqIO.parse(file_path, 'genbank'):
            for feature in rec.features:
                if feature.type == 'CDS':
                    if feature.qualifiers.get('gene') and feature.qualifiers.get('translation'):
                        if feature.qualifiers['gene'][0].lower() not in b:
                            b.append(feature.qualifiers['gene'][0].lower())
                    else:
                        print(f"{feature} in genbank file has no gene name or translation!, please verify!!!")
        # lower adj
        for x in range(len(a) - 1, -1, -1):
            if a[x].lower() not in [y.lower() for y in b]:
                a.remove(a[x])
    # save gene to gene results
    saved_gene = os.path.join(os.path.dirname(work_dir), 'pro_gene.txt')
    gene_name_file = open((saved_gene), 'w')
    for gene_name in a:
        # trans gene name to rec name
        for rec in SeqIO.parse(ref_file, format='genbank'):
            for feature in rec.features:
                if feature.type == 'CDS':
                    if feature.qualifiers['gene'][0].lower() == gene_name:
                        gene_name_file.write(f"{feature.qualifiers['gene'][0]}\n")
                        break

    print(f"Total {len(a)} common protein genes, which have been saved in: \n\t\t{os.path.abspath(saved_gene)}\n"
          f"{'-'*80}")
    gene_name_file.close()
    # extract cds to each fasta file
    fasta_file = os.path.join(output_dir, 'pro_fasta')
    if os.path.exists(fasta_file):
        print(f"The file path has existed, please change the output directory:\n\t\t{os.path.abspath(fasta_file)}\n"
              f"{'-'*80}")
    else:
        os.makedirs(fasta_file)
        for j in os.listdir(work_dir):
            if j.endswith('gb') or j.endswith('gbk'):
                gb_files = os.path.join(work_dir, j)
                fasta_files = os.path.join(fasta_file, str(j.split('.')[0]) + '.fasta')
                with open(fasta_files, 'w') as gg:
                    gg.write(f">{j.split('.')[0]}\n")
                    for x in a:
                        for rec in SeqIO.parse(gb_files, format="genbank"):
                            seqs = []
                            for feature in rec.features:
                                if feature.type == "CDS":
                                    if feature.qualifiers['gene'][0].lower() == x.lower():
                                        seqs.append(feature.qualifiers['translation'][0])
                            if len(seqs) == 1:
                                gg.write(str(seqs[0]))
                            if len(seqs) == 2:
                                seqs.remove(seqs[0]) if len(seqs[0]) <= len(seqs[1]) else seqs.remove(seqs[1])
                                gg.write(str(seqs[0]))
                gg.close()
    # merge sequences to one fasta_file
    merge_file = os.path.join(output_dir, 'merge_pro.fasta')
    merge_fasta = open(merge_file, 'w')
    for fa_list in os.listdir(fasta_file):
        fa_file = os.path.join(fasta_file, fa_list)
        for rec in SeqIO.parse(fa_file, format='fasta'):
            merge_fasta.write(f">{rec.id}\n{rec.seq}\n")
    merge_fasta.close()
    print(f"The merged protein file is saved in:\n\t\t {os.path.abspath(merge_file)}\n"
          f"{'-'*80}")


def cds_for_phy(ref_file, output_dir):
    a = []
    work_dir = os.path.dirname(ref_file)
    for rec in SeqIO.parse(ref_file, format='genbank'):
        for feature in rec.features:
            if feature.type == 'CDS':
                if feature.qualifiers.get('gene'):
                    if feature.qualifiers['gene'][0].lower() not in a:
                        a.append(feature.qualifiers['gene'][0].lower())
                else:
                    print(f"{feature} in genbank file has no gene name!, please verify!!!")
    # delete not common cds
    for i in os.listdir(work_dir):
        file_path = os.path.join(work_dir, i)
        b = []
        for rec in SeqIO.parse(file_path, 'genbank'):
            for feature in rec.features:
                if feature.type == 'CDS':
                    if feature.qualifiers.get('gene'):
                        if feature.qualifiers['gene'][0].lower() not in b:
                            b.append(feature.qualifiers['gene'][0].lower())
                    else:
                        print(f"{feature} in genbank file has no gene name!, please verify!!!")
        # lower adj
        for x in range(len(a) - 1, -1, -1):
            if a[x].lower() not in [y.lower() for y in b]:
                a.remove(a[x])
    # save gene to gene results
    saved_gene = os.path.join(os.path.dirname(work_dir), 'cds_gene.txt')
    gene_name_file = open((saved_gene), 'w')
    for gene_name in a:
        # trans gene name to rec name
        for rec in SeqIO.parse(ref_file, format='genbank'):
            for feature in rec.features:
                if feature.type == 'CDS':
                    if feature.qualifiers['gene'][0].lower() == gene_name:
                        gene_name_file.write(f"{feature.qualifiers['gene'][0]}\n")
                        break

    print(f"Total {len(a)} common cds genes, which have been saved in: \n\t\t{os.path.abspath(saved_gene)}\n"
          f"{'-'*80}")
    gene_name_file.close()
    # extract cds to each fasta file
    fasta_file = os.path.join(output_dir, 'cds_fasta')
    if os.path.exists(fasta_file):
        print(f"The file path has existed, please change the output directory:\n\t\t{os.path.abspath(fasta_file)}\n"
              f"{'-'*80}")
    else:
        os.makedirs(fasta_file)
        for j in os.listdir(work_dir):
            if j.endswith('gb') or j.endswith('gbk'):
                gb_files = os.path.join(work_dir, j)
                fasta_files = os.path.join(fasta_file, str(j.split('.')[0]) + '.fasta')
                with open(fasta_files, 'w') as gg:
                    gg.write(f">{j.split('.')[0]}\n")
                    for x in a:
                        for rec in SeqIO.parse(gb_files, format="genbank"):
                            seqs = []
                            for feature in rec.features:
                                if feature.type == "CDS":
                                    if feature.qualifiers['gene'][0].lower() == x.lower():
                                        seqs.append(feature.extract(rec.seq))
                            if len(seqs) == 1:
                                gg.write(str(seqs[0]))
                            if len(seqs) == 2:
                                seqs.remove(seqs[0]) if len(seqs[0]) <= len(seqs[1]) else seqs.remove(seqs[1])
                                gg.write(str(seqs[0]))
                gg.close()
    # merge sequences to one fasta_file
    merge_file = os.path.join(output_dir, 'merge_cds.fasta')
    merge_fasta = open(merge_file, 'w')
    for fa_list in os.listdir(fasta_file):
        fa_file = os.path.join(fasta_file, fa_list)
        for rec in SeqIO.parse(fa_file, format='fasta'):
            merge_fasta.write(f">{rec.id}\n{rec.seq}\n")
    merge_fasta.close()
    print(f"The merged cds file is saved in:\n\t\t {os.path.abspath(merge_file)}\n"
          f"{'-'*80}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-d", "--input_dir", help="Input the directory of genbank files")
    parser.add_argument("-m", "--mode", choices=['cds', 'pro'],
                        help="Mode: cds for common cds sequences; pro for common protein sequences")
    args = parser.parse_args()
    if args.input_dir:
        ref_file = find_ref(args.input_dir)
        if ref_file and (ref_file.endswith('gb') or ref_file.endswith('gbk')):
            file_path = os.path.abspath(ref_file)
            if args.mode == 'cds':
                save_path = os.path.join(os.path.dirname(args.input_dir), 'cds')
                cds_for_phy(file_path, save_path)
            elif args.mode == 'pro':
                save_path = os.path.join(os.path.dirname(args.input_dir), 'pro')
                pro_for_phy(file_path, save_path)
            else:
                print("Two modes are provided, you must specify one of them")
        else:
            print(f"{ref_file} is not genbank format files, Please ends with 'gb' or 'gbk'\n{'-' * 40}\n")
    else:
        parser.print_help()
