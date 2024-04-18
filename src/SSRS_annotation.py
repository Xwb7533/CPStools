import os
from Bio import SeqIO
import re
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Find SSRs in chloroplast genomes.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com
        """

        return f"{description}\n\n{help_text}"


def check_files(input_file):
    try:
        if input_file.endswith('.fasta') or input_file.endswith('.fa'):
            for rec in SeqIO.parse(input_file, 'fasta'):
                print("Please input Genbank format file")
        if input_file.endswith('.gb') or input_file.endswith('.gbk'):
            for rec in SeqIO.parse(input_file, 'genbank'):
                return str(rec.seq)
        raise ValueError(f"Unsupported file format for '{input_file}'. "
                         f"Please use FASTA (.fasta or .fa) or GenBank (.gb or .gbk) files.")
    except Exception as e:
        if isinstance(e, FileNotFoundError):
            raise  # Reraise the FileNotFoundError
        elif os.path.exists(input_file):
            raise ValueError(f"Error processing file '{input_file}': {e}")
        else:
            raise FileNotFoundError(f"No such file: {input_file}")


def find_SSRs(input_file, *length):
    type_length = []
    if len(length) == 0:
        type_length = [10, 6, 5, 4, 4, 4]
    elif len(length) == 6:
        type_length = length
    else:
        print("Please input type length equal 6 or 0")
        return []
    k1, k2, k3, k4, k5, k6 = type_length

    print(f"The parameter is set as:\n"
          f"Mononucleotide:{k1}\n"
          f"Dinucleotide:{k2}\n"
          f"Trinucleotide:{k3}\n"
          f"Tetranucleotide:{k4}\n"
          f"Pentanucleotide:{k5}\n"
          f"Hexanucleotide:{k6}\n")
    my_seq = check_files(input_file)
    matches = {
        'match1': re.finditer(r'([ATCG])\1{{{0},}}'.format(k1 - 1), my_seq),
        'match2': re.finditer(r'(((?!AA|TT|CC|GG)[ATCG]){{2}})\1{{{0},}}'.format(k2 - 1), my_seq),
        'match3': re.finditer(r'(((?!AAA|TTT|CCC|GGG)[ATCG]){{3}})\1{{{0},}}'.format(k3 - 1), my_seq),
        'match4': re.finditer(r'(?!(\w)(\w)\1\2)(((?!AAAA|TTTT|CCCC|GGGG)[ATCG]){{4}})\3{{{0},}}'.format(k4 - 1),
                              my_seq),
        'match5': re.finditer(r'(((?!AAAAA|TTTTT|CCCCC|GGGGG)[ATCG]){{5}})\1{{{0},}}'.format(k5 - 1), my_seq),
        'match6': re.finditer(
            r'(?!(\w)(\w)(\w)\1\2\3)(((?!AAAAAA|TTTTTT|CCCCCC|GGGGGG)[ATCG]){{6}})\4{{{0},}}'.format(k6 - 1), my_seq)
    }
    all_matches = []
    for index_ in range(1, 7):
        group_name = 'match' + str(index_)
        for m in matches[group_name]:
            start = m.start() + 1
            if not any(start == int(match.split('\t')[2]) for match in all_matches):
                character = m.group()[:index_]
                length = len(m.group()) // index_
                end = m.end()
                output = f"{character}\t{length}\t{start}\t{end}"

                # print(output)
                all_matches.append(output)

    all_matches_sorted = sorted(all_matches, key=lambda x: int(x.split('\t')[2]))
    file_name = os.path.basename(input_file).split('.')[0]
    save_name = file_name + "_SSRs_results.txt"
    # save_name2 = file_name + "_SSRs_simplification.txt"
    file_save = os.path.join(os.path.dirname(input_file), save_name)
    with open(file_save, 'w') as ff:
        ff.write("type\tlength\tstart\tend\n")
        for SSRs_ in all_matches_sorted:
            ff.write(f"{SSRs_}\n")
    print(f"results:\ttotal {len(all_matches)} SSRs were detected!\n{'-' * 80}")


def IGS_extract(input_file):
    for rec in SeqIO.parse(input_file, format='genbank'):
        genome_length = [[int(part.end)] for part in rec.features[0].location.parts][0][0]
        my_seq = rec.seq
        all_feature = []
        all_info = []
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                all_feature.append(feature)
        for i in range(len(all_feature)):
            gene_name = all_feature[i].qualifiers['gene'][0]
            gene_location = all_feature[i].location.parts
            gene1_exon_info = [[int(part.start), int(part.end), part.strand] for part in gene_location]
            exon_number = len(gene_location)
            if exon_number == 1:
                all_info.append(f"{gene_name}\t{gene1_exon_info[0][0]}\t{gene1_exon_info[0][1]}\t{gene1_exon_info[0][2]}")
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
        save_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + '_SSRs.txt')
        save_file_w = open(save_file, 'w')
        for i in range(len(all_info)-1):
            info_list = all_info[i].split('\t')
            next_list = all_info[i+1].split('\t')
            #print(info_list, next_list, info_list[0].split('_')[0], next_list[0])
            save_file_w.write(f"{info_list[0]}\t{info_list[1]}\t{info_list[2]}\tGene\n")
            save_file_w.write(f"{info_list[0]}-{next_list[0]}\t{info_list[-2]}\t{next_list[1]}\tIGS\n")
        save_file_w.write(f"{all_info[-1][:-1]}\tGene\n")
        end_gene_info = all_info[-1].split('\t')
        start_gene_info = all_info[0].split('\t')
        if int(end_gene_info[-2]) < int(start_gene_info[1]):
            save_file_w.write(f"{end_gene_info[0]}-{start_gene_info}[0]\t{end_gene_info[-2]}\t{start_gene_info[1]}"
                              f"\tGene\n")
        else:
            if int(end_gene_info[2]) < genome_length:
                save_file_w.write(f"{end_gene_info[0]}-{start_gene_info[0]}\t{end_gene_info[-2]}\t{genome_length}\t0\t\
                {start_gene_info[1]}\tIGS\n")
            else:
                pass
        save_file_w.close()

        # change location type 
        output_ = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + '_SSRs2.txt')
        change_file = open(save_file, 'r').readlines()
        with open(output_, 'w') as gg:
            for loc2 in change_file:
                loc_type = loc2.split('\t')[0]
                if "matK" in loc_type and "trnK-UUU" in loc_type:
                    loc2_list = loc2.split('\t')
                    loc2_list[-1] = "Intron"
                    gg.write('\t'.join(loc_x for loc_x in loc2_list) + '\n')
                else:
                    if len(loc_type.split('_')) == 3:
                        N1, N3 = loc2.split('_')[0], loc2.split('_')[1]
                        if N1 in N3:
                            loc2_list = loc2.split('\t')
                            loc2_list[-1] = "Intron"
                            gg.write('\t'.join(loc_x for loc_x in loc2_list) + '\n')
                        else:
                            gg.write(loc2)
                    else:
                        gg.write(loc2)


        # get location
        file_name = os.path.basename(input_file).split('.')[0]
        save_name = file_name + "_SSRs_results.txt"
        file1 = os.path.join(os.path.dirname(input_file), save_name)
        st_lines = open(output_, 'r').readlines()[1:]

        final_name = file_name + "_SSRs_loc_results.txt"
        final_results = os.path.join(os.path.dirname(input_file), final_name)
        data = []
        st_lines2 = open(file1, 'r').readlines()[1:]
        for line in st_lines:
            parts = line.strip().split('\t')
            name, start, end, loc = parts[0], int(parts[1]), int(parts[2]), parts[-1]
            data.append((name, start, end, loc))

        printed_combinations = set()
        # write output results
        with open(final_results, 'w') as ff:
            ff.write("type\tlength\tstart\tend\tloc\tloc_type\n")
            for i in st_lines2:
                st_list = i.split('\t')
                st_join = '\t'.join(info.strip() for info in st_list)
                target_list = list(range(int(st_list[2]), int(st_list[3])))
                for num in target_list:
                    for name, start, end, loc in data:
                        if start <= num <= end:
                            combination = (tuple(target_list), name)
                            if combination not in printed_combinations:
                                ff.write(f"{st_join}\t{name}\t{loc}\n")
                                printed_combinations.add(combination)
                            break
    os.remove(output_)
    os.remove(save_file)
    os.remove(file1)
    print(f"The results was written into:\n\t\t{os.path.abspath(final_results)}\n{'-' * 80}")
                        


def main():
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--input_file', help='GenBank format file', required=True)
    parser.add_argument('-k', '--kmer_length', help='SSRs length, default is 10,6,5,4,4,4')
    # parser.add_argument('-o', '--output', help='output file')

    args = parser.parse_args()
    if args.input_file:
        try:
            file_path = os.path.abspath(args.input_file)
            kmer_lengths = list(map(int, args.kmer_length.split(','))) if args.kmer_length else []
            check_results = check_files(file_path)
            if check_results:
                find_SSRs(file_path, *kmer_lengths)
                IGS_extract(file_path)
        except (ValueError, FileNotFoundError) as e:
            print(e)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

