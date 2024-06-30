import os
from Bio import SeqIO
import argparse
import sys
import subprocess


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Calculate Pi valus from Genbank files and sort as cp order.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""

        return f"{description}\n\n{help_text}"


def IGS_extract(input_file, fasta_dir, info_dir):
    for rec in SeqIO.parse(input_file, format='genbank'):
        genome_length = [[int(part.end)] for part in rec.features[0].location.parts][0][0]
        my_seq = rec.seq
        all_feature = []
        all_info = []
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                all_feature.append(feature)
        for i in range(len(all_feature)):
            #print(input_file)
            # print(all_feature[i])
            gene_name = all_feature[i].qualifiers['gene'][0]
            gene_location = all_feature[i].location.parts
            gene1_exon_info = [[int(part.start), int(part.end), part.strand] for part in gene_location]
            exon_number = len(gene_location)
            if exon_number == 1:
                all_info.append(f"{gene_name}\t{gene1_exon_info[0][0]}\t{gene1_exon_info[0][1]}\t{gene1_exon_info[0][2]}")
            if exon_number == 2:
                if gene1_exon_info[0][1] == genome_length:
                    all_info.append(f"{gene_name}\t{gene1_exon_info[0][0]}\t{gene1_exon_info[0][1]}\t{gene1_exon_info[1][0]}\t{gene1_exon_info[1][1]}\t{gene1_exon_info[0][2]}")
                else:
                    all_info.append(f"{gene_name}_1\t{gene1_exon_info[0][0]}\t{gene1_exon_info[0][1]}\t{gene1_exon_info[0][2]}")
                    all_info.append(f"{gene_name}_2\t{gene1_exon_info[1][0]}\t{gene1_exon_info[1][1]}\t{gene1_exon_info[1][2]}")
            if exon_number == 3:
                all_info.append(f"{gene_name}_1\t{gene1_exon_info[0][0]}\t{gene1_exon_info[0][1]}\t{gene1_exon_info[0][2]}")
                all_info.append(f"{gene_name}_2\t{gene1_exon_info[1][0]}\t{gene1_exon_info[1][1]}\t{gene1_exon_info[1][2]}")
                all_info.append(f"{gene_name}_3\t{gene1_exon_info[2][0]}\t{gene1_exon_info[2][1]}\t{gene1_exon_info[2][2]}")
        all_info = list(set(all_info))
        all_info.sort(key=lambda x: int(x.split('\t')[1]))
        save_file = os.path.join(info_dir, os.path.basename(input_file).split('.')[0] + '_intergenic_location.txt')
        save_file_w = open(save_file, 'w')
        for i in range(len(all_info)-1):
            info_list = all_info[i].split('\t')
            next_list = all_info[i+1].split('\t')
            save_file_w.write(f"{info_list[0]}-{next_list[0]}\t{info_list[-2]}\t{next_list[1]}\n")
        end_gene_info = all_info[-1].split('\t')
        start_gene_info = all_info[0].split('\t')
        if int(end_gene_info[-2]) < int(start_gene_info[1]):
            save_file_w.write(f"{end_gene_info[0]}-{start_gene_info}[0]\t{end_gene_info[-2]}\t{start_gene_info[1]}\n")
        else:
            if int(end_gene_info[2]) < genome_length:
                save_file_w.write(f"{end_gene_info[0]}-{start_gene_info[0]}\t{end_gene_info[-2]}\t{genome_length}\t0\t\
                {start_gene_info[1]}\n")
            else:
                pass
        save_file_w.close()
        all_fasta_file = os.path.join(fasta_dir, os.path.basename(input_file).split('.')[0] + '_IGS.fasta')
        all_fasta = open(all_fasta_file, 'w')
        save_results = open(save_file, 'r')
        result_line = save_results.readline().strip()
        while result_line:
            result_line_list = result_line.split('\t')
            if len(result_line_list) == 3:
                if int(result_line_list[2]) > int(result_line_list[1]):
                    all_fasta.write(f">{result_line_list[0]}\n{my_seq[int(result_line_list[1]):int(result_line_list[2])]}\n")
                    result_line = save_results.readline().strip()
                else:
                    # print(f"{result_line_list[0]} has overlap!")
                    result_line = save_results.readline().strip()

            else:
                all_fasta.write(f">{result_line_list[0]}\n{my_seq[int(result_line_list[1]):int(result_line_list[2])]}\
                {my_seq[int(result_line_list[3]):int(result_line_list[4])]}\n")
                result_line = save_results.readline().strip()
        all_fasta.close()


def common_IGS(input_dir):
    all_common = []
    input_file = None
    
    # Find the first .fasta file in the input_dir
    for file in os.listdir(input_dir):
        if file.endswith('.fasta'):
            input_file = os.path.join(input_dir, file)
            break
    
    if input_file is None:
        print("No fasta file found in the directory.")
        return
    
    # Parse the input file to get the initial list of sequence IDs
    for rec in SeqIO.parse(input_file, format='fasta'):
        all_common.append(rec.id)
    
    work_dir = os.path.dirname(input_file)
    
    # Process each .fasta file in the directory to find common sequence IDs
    for fasta_file in os.listdir(work_dir):
        if fasta_file.endswith('.fasta'):
            single_IGS = []
            fasta_path = os.path.join(work_dir, fasta_file)
            
            for rec2 in SeqIO.parse(fasta_path, format='fasta'):
                single_IGS.append(rec2.id)
            
            all_common = [common_id for common_id in all_common if common_id.lower() in (single_id.lower() for single_id in single_IGS)]
    
    # Save the common sequence IDs to a file
    with open(os.path.join(work_dir, 'cp_sort_IGS.txt'), 'w') as cp_sort_IGS_file:
        for i in all_common:
            cp_sort_IGS_file.write(f"{i}\n")
    
    
    save_dir = os.path.join(work_dir, 'unalign_common_IGS')
    print(f"The common intergenic fasta number is {len(all_common)}, and the unalignment files are saved at \n\t\t\t{save_dir}")
    os.mkdir(save_dir)
    
    # Save the common sequences to separate files
    for common_name in all_common:
        save_file_path = os.path.join(save_dir, common_name + '.fasta')
        
        with open(save_file_path, 'w') as save_file:
            for fasta_file in os.listdir(work_dir):
                if fasta_file.endswith('.fasta'):
                    fasta_path = os.path.join(work_dir, fasta_file)
                    
                    for rec in SeqIO.parse(fasta_path, format='fasta'):
                        if rec.id.lower() == common_name.lower():
                            save_file.write(f">{fasta_file.split('_IGS')[0]}\n{rec.seq}\n")
    
    return save_dir


def common_gene_extract(input_dir):
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"The input directory {input_dir} does not exist.")
    
    # 获取目录中的第一个文件作为参考
    reference_file = None
    for file in os.listdir(input_dir):
        if file.endswith('.gb') or file.endswith('.gbk'):
            reference_file = os.path.join(input_dir, file)
            break
    
    if reference_file is None:
        raise FileNotFoundError("No GenBank files found in the input directory.")
    
    work_dir = input_dir
    all_gene = []
    
    # 使用参考文件提取基因列表
    for rec in SeqIO.parse(reference_file, format='genbank'):
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                if feature.qualifiers['gene'][0] not in all_gene:
                    all_gene.append(feature.qualifiers['gene'][0])
    
    # 遍历目录中的所有文件，提取公共基因
    for files in os.listdir(work_dir):
        if files.endswith('gb') or files.endswith('gbk'):
            single_gene = []
            gb_file = os.path.join(work_dir, files)
            for rec in SeqIO.parse(gb_file, format='genbank'):
                for feature in rec.features:
                    if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                        if feature.qualifiers['gene'][0] not in single_gene:
                            single_gene.append(feature.qualifiers['gene'][0])
            # 删除不公共的基因
            for gene_index in range(len(all_gene)-1, -1, -1):
                if all_gene[gene_index].lower() not in [y.lower() for y in single_gene]:
                    all_gene.remove(all_gene[gene_index])
    
    gene_name_file = os.path.join(os.path.dirname(work_dir), 'gene_cp_sort.txt')
    # 将结果保存到文件中
    with open(gene_name_file, 'w') as f:
        for gene in all_gene:
            f.write(f"{gene}\n")
    # return gene_name_file
    
    save_dir = os.path.join(os.path.dirname(work_dir), 'common_gene')
    if os.path.exists(save_dir):
        print(
            f"\t\t########  Run failed !!!  ########\n"
            f"The save directory has been existed, please delete the directory!\n"
            f"\t\t\t{os.path.abspath(save_dir)}")
        sys.exit()
    os.mkdir(save_dir)
    print(f"The common gene fasta number is {len(all_gene)}, and the unalignment files are saved at \
        \n\t\t\t{save_dir}")
    for gene_name in all_gene:
        file_name = str(gene_name) + '.fasta'
        file_path = os.path.join(save_dir, file_name)
        with open(file_path, 'w') as fasta_file:
            for gb_file in os.listdir(work_dir):
                if gb_file.endswith('gb') or gb_file.endswith('gbk'):
                    gb_file_path = os.path.join(work_dir, gb_file)
                    fasta_file.write(f">{gb_file.split('.')[0]}\n")
                    for rec in SeqIO.parse(gb_file_path, format='genbank'):
                        my_seqs = []
                        for feature in rec.features:
                            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                                if feature.qualifiers['gene'][0].lower() == gene_name.lower():
                                    my_seqs.append(feature.extract(rec.seq))
                        if len(my_seqs) == 1:
                            fasta_file.write(f"{my_seqs[0]}\n")
                        if len(my_seqs) == 2:
                            my_seqs.remove(my_seqs[0]) if len(my_seqs[0]) <= len(my_seqs[1]) else my_seqs.remove(my_seqs[1])
                            fasta_file.write(f"{my_seqs[0]}\n")
    return save_dir



def is_mafft_available():
    """Check if MAFFT is available in the system path."""
    try:
        subprocess.run(["mafft", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def multi_mafft(input_dir, mafft_path=None):
    if not is_mafft_available() and mafft_path is None:
        raise ValueError("MAFFT is not in the environment path and no MAFFT path is provided.")
    
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"The input directory {input_dir} does not exist.")
    
    align_dir = os.path.join(os.path.dirname(input_dir), 'align_gene')
    os.makedirs(align_dir, exist_ok=True)

    mafft_executable = mafft_path if mafft_path else 'mafft'
    
    for file in os.listdir(input_dir):
        input_file = os.path.join(input_dir, file)
        output_file = os.path.join(align_dir, f'{file}')
        
        command = f'{mafft_executable} --auto {input_file} > {output_file}'
        with open(os.devnull, 'w') as devnull:
            result = subprocess.run(command, shell=True, stdout=devnull, stderr=devnull)
    return align_dir

def calculate_Pi_values(work_dir):
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
            # print(f"{align_fasta[:-6]}\t{final_pi}")
    
    Pi_results = os.path.join(work_dir, 'Pi_results.txt')
    with open(Pi_results, 'w') as ff:
        for each_pi in all_pi_results:
            ff.write(f"{each_pi}\n")
    return Pi_results

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
                # print(IGS1, end='')
                reuslts_file.write(IGS1)
    reuslts_file.close()
    return results_file_path

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
                # print(IGS1, end='')
                reuslts_file.write(IGS1)
    reuslts_file.close()
    return results_file_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-d", "--work_dir", required=True, help="Input directory of genbank files")
    parser.add_argument("-m", "--mafft_path", help="Path to MAFFT executable if not in environment path")
    args = parser.parse_args()
    if not args.work_dir:
        parser.error("Please provide the input directory of genbank files using -i or --work_dir")
    else:
        save_results_dir = common_gene_extract(os.path.abspath(args.work_dir))
        # print(save_results_dir)
        # extract fasta
        fasta_dir = os.path.abspath(os.path.join(os.path.dirname(args.work_dir), 'IGS'))
        info_dir = os.path.join(os.path.dirname(args.work_dir), 'info')
        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(info_dir, exist_ok=True)
        file_name = ''
        for i in os.listdir(args.work_dir):
            if i.endswith('gb') or i.endswith('gbk'):
                file_path = os.path.join(args.work_dir, i)
                file_name = os.path.join(fasta_dir, os.path.basename(i).split('.')[0] + '_IGS.fasta')
                # print(file_path)
                IGS_extract(file_path, fasta_dir, info_dir)
            else:
                print(f"please input genbank format files and endwith 'gb' or 'gbk'! ")
        
        # print(f"file name is {file_name}")
        if is_mafft_available() or args.mafft_path:
            IGS_results_dir = multi_mafft(common_IGS(fasta_dir), args.mafft_path)
            print('-'*80)
            print("Performing multiple sequence alignment, please wait......")
            gene_results_dir = multi_mafft(save_results_dir, args.mafft_path)
            print('-'*80)
            print(f'The aligned IGS-sequences are saved at {IGS_results_dir}')
            print('-'*80)
            print(f'The aligned IGS-sequences are saved at {gene_results_dir}')

            # # 计算 Pi 值
            IGS_Pi_values = calculate_Pi_values(IGS_results_dir)
            gene_Pi_values = calculate_Pi_values(gene_results_dir)
            gene_name_file = os.path.abspath(os.path.join(os.path.dirname(args.work_dir), 'gene_cp_sort.txt'))
            IGS_name_file = os.path.abspath(os.path.join(os.path.dirname(args.work_dir), 'IGS/cp_sort_IGS.txt'))
            # print(IGS_Pi_values, gene_Pi_values, gene_name_file)
            fianl_IGS = IGS_sort_as_cp_order(IGS_Pi_values, IGS_name_file)
            final_gene = gene_sort_as_cp_order(gene_Pi_values, gene_name_file)
            print('-'*80)
            print(f"Pi values of IGS sorted as cp order are saved in: {fianl_IGS}")
            print('-'*80)
            print(f"Pi values of gene sorted as cp order are saved in: {final_gene}")
        else:
            print("Please provide the abspath of mafft")
            

