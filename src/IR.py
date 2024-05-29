import os, sys, re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """Find four regions in chloroplast genomes.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""

        return f"{description}\t\t{help_text}"


def check_files(input_file):
    try:
        if input_file.endswith('.fasta') or input_file.endswith('.fa'):
            for rec in SeqIO.parse(input_file, 'fasta'):
                if not re.fullmatch(r'[ATCG]*', str(rec.seq)):
                    raise ValueError("Sequence contains bases other than A, T, C, G.")
                    sys.exit()
                else:
                    return rec.seq
        if input_file.endswith('.gb') or input_file.endswith('.gbk'):
            for rec in SeqIO.parse(input_file, 'genbank'):
                if not re.fullmatch(r'[ATCG]*', str(rec.seq)):
                    raise ValueError("Sequence contains bases other than A, T, C, G.")
                    sys.exit()
                else:
                    return rec.seq
        raise ValueError(f"Unsupported file format for '{input_file}'. "
                         f"Please use FASTA (.fasta or .fa) or GenBank (.gb or .gbk) files.")
    except Exception as e:
        if isinstance(e, FileNotFoundError):
            raise  # Reraise the FileNotFoundError
        elif os.path.exists(input_file):
            raise ValueError(f"Error processing file '{input_file}': {e}")
        else:
            raise FileNotFoundError(f"No such file: {input_file}")


def find_repeat_regions(sequences):
    replace_dit = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C",
        "a": "t",
        "c": "g",
        "t": "a",
        "g": "c"
    }
    sequences = str(Seq(str(sequences))).lower()
    rev_seq = str(Seq(sequences).reverse_complement())
    # seed test
    all_list = []
    for i in range(0, len(sequences), 500):
        if i + 500 <= len(sequences):
            test_seq = sequences[i:i + 500]
            if test_seq in rev_seq:
                all_list.append([i, i + 500])
    # print(all_list)
    if not all_list:
        print("No repeated sequenced longer than 1,000 bp was detected!")
    else:
        merged_sublists = []
        start = all_list[0][0]
        end = all_list[0][1]

        for sublist in all_list[1:]:
            sublist_start = sublist[0]
            sublist_end = sublist[1]

            if sublist_start == end:
                end = sublist_end
            else:
                merged_sublists.append([start, end])
                start = sublist_start
                end = sublist_end

        merged_sublists.append([start, end])
        # print(merged_sublists)
        # 找出最长的重复序列
        longest_sublist = max(merged_sublists, key=lambda x: x[1] - x[0])
        # print(f"longest_sublist: {longest_sublist}")

        # 找到最长序列的反向互补序列的位置
        start_pos, start_end = longest_sublist[0], longest_sublist[1]
        max_seq = str(Seq(sequences[start_pos:start_end]).reverse_complement())
        start_loc = sequences.find(max_seq)
        if start_loc != -1:
            seq_end = start_loc + len(max_seq)
            start_list = [start_loc, start_loc + len(max_seq)]
        else:
            print("The assembled may be wrong!")
            sys.exit()

        # 判断前后关系

        def sort_lists(list1, list2):
            if list1[0] < list2[0]:
                return [list1, list2]
            else:
                return [list2, list1]

        _list = sort_lists(start_list, longest_sublist)
        sorted_list = [item for sublist in _list for item in sublist]
        # print(sorted_list)

        # 延伸并确定出重复区域的序列

        start1, end1, start2, end2 = sorted_list[0], sorted_list[1], sorted_list[2], sorted_list[3]
        # print(start1, end1, start2, end2)

        # 内侧延伸
        while sequences[end1 - 1] == replace_dit[sequences[start2]]:
            end1 += 1
            start2 -= 1

        # print(start1, end1, start2, end2)
        # IR 不在1处
        if start1 != 0:
            if end2 != len(sequences):
                # 继续延伸
                if sequences[start1 - 1] == replace_dit[sequences[end2]]:
                    new_start = start1 - 1
                    new_end = end2
                    while new_end <= len(sequences) - 1 and sequences[new_start] == replace_dit[sequences[new_end]]:
                        new_start -= 1
                        new_end += 1
                    # print("******")
                    # print(new_start, new_end)
                    if new_end == len(sequences):
                        # 跨区域延伸
                        if sequences[new_start] == replace_dit[sequences[0]]:
                            new_end = 1
                            new_start -= 1
                            while sequences[new_start] == replace_dit[sequences[new_end]]:
                                new_start -= 1
                                new_end += 1
                            # LSC在前
                            if (start2 - end1) <= (new_start - new_end):
                                result = (
                                    f"LSC:{new_end + 1}-{new_start + 1}\tIRb:{new_start + 2}-{end1 - 1}\t"
                                    f"SSC:{end1}-{start2 + 1}\tIRa:{start2 + 2}-{len(sequences)},1-{new_end}"
                                )
                            # SSC在前
                            else:
                                result = (
                                    f"LSC:{end1}-{start2 + 1}\tIRb:{start2 + 2}-{len(sequences)},1-{new_end}\t"
                                    f"SSC:{new_end + 1}-{new_start + 1}\tIRa:{start1 + 2}-{end1 - 1}"
                                )
                        # 不跨区域延伸
                        else:
                            # LSC在前
                            if (start2 - end1) <= new_start:
                                result = (
                                    f"LSC:{1}-{new_start + 1}\tIRb:{new_start + 2}-{end1 - 1}\t"
                                    f"SSC:{end1}-{start2 + 1}\tIRa:{start2 + 2}-{new_end}"
                                )
                            # SSC在前
                            else:
                                result = (
                                    f"LSC:{end1}-{start2 + 1}\tIRb:{start2 + 2}-{new_end}\t"
                                    f"SSC:1-{new_start + 1}\tIRa:{new_start + 2}-{end1 - 1}"
                                )
                    # 不跨区域
                    else:
                        # LSC在前
                        if (start2 - end1) <= (len(sequences) - new_end + new_start):
                            result = (
                                f"LSC:{new_end + 1}-{len(sequences)},1-{new_start + 1}\tIRb:{new_start + 2}-{end1 - 1}\t"
                                f"SSC:{end1}-{start2 + 1}\tIRa:{start2 + 2}-{new_end}"
                            )
                        # SSC在前
                        else:
                            result = (
                                f"LSC:{end1}-{start2 + 1}\tIRb:{start2 + 2}-{new_end}\t"
                                f"SSC:{new_end + 1}-{len(sequences)},1-{new_start + 1}\tIRa:{start1 + 2}-{end1 - 1}"
                            )
                # 停止延伸
                else:
                    # LSC在前
                    if (start2 - end1) <= (len(sequences) - end2 + start1):
                        result = (
                            f"LSC:{end2 + 1}-{len(sequences)},1-{start1 + 1}\tIRb:{start1 + 2}-{end1 - 1}\t"
                            f"SSC:{end1}-{start2 + 1}\tIRa:{start2 + 2}-{end2}"
                        )
                    # SSC在前
                    else:
                        result = (
                            f"LSC:{end1}-{start2 + 1}\tIRb:{start2 + 2}-{end2}\t"
                            f"SSC:{end2 + 1}-{len(sequences)},1-{start1 + 1}\tIRa:{start1 + 2}-{end1 - 1}"
                        )
            # 直接跨区域
            else:
                # 延伸
                if sequences[start1 - 1] == replace_dit[sequences[0]]:
                    new_end = 0
                    new_start = start1 - 1
                    while sequences[new_start] == replace_dit[sequences[new_end]]:
                        new_start -= 1
                        new_end += 1
                    # LSC在前
                    if (start2 - end1) <= (new_start - new_end):

                        result = (
                            f"LSC:{new_end + 1}-{new_start + 1}\tIRb:{new_start + 2}-{end1 - 1}\t"
                            f"SSC:{end1}-{start2 + 1}\tIRa:{start2 + 2}-{len(sequences)},1-{new_end}"
                        )
                    # SSC在前
                    else:
                        result = (
                            f"LSC:{end1}-{start2 + 1}\tIRb:{start2 + 2}-{len(sequences)},1-{new_end}\t"
                            f"SSC:{new_end + 1}-{new_start}\tIRa:{start1 + 1}-{end1 - 1}"
                        )
                # 不延伸
                else:
                    if (start2 - end1) <= start1:
                        result = (
                            f"LSC:1-{start1}\tIRb:{start1 + 1}-{end1 - 1}\t"
                            f"SSC:{end1}-{start2 + 1}\tIRa:{start2 + 2}-{end2}"
                        )
                    else:
                        result = (
                            f"LSC:{end1}-{start2 + 1}\tIRb:{start2 + 2}-{end2}\t"
                            f"SSC:{end2 + 1}-{len(sequences)},1-{start1}\tIRa:{start1 + 1}-{end1 - 1}"
                        )

        # IR 区域在1处
        else:
            # 跨区域延伸
            if sequences[-1] == replace_dit[sequences[end2]]:
                new_start = len(sequences) - 1
                new_end = end2
                while sequences[new_start] == replace_dit[sequences[new_end]]:
                    new_start -= 1
                    new_end += 1
                # 判断LSC 和 SSC
                # LSC在前
                if (start2 - end1) <= (new_start - new_end):
                    result = (
                        f"LSC:{new_end + 1}-{new_start + 1}\tIRb:{new_start + 2}-{len(sequences)},1-{end1 - 1}\t"
                        f"SSC:{end1}-{start2 + 1}\tIRa:{start2 + 2}-{end2}"
                    )
                # SSC在前
                else:
                    result = (
                        f"LSC:{end1}-{start2 + 1}\tIRb:{start2 + 2}-{new_end}\tSSC{new_end + 1}-{new_start + 1}\t"
                        f"IRa:{new_start + 2}-{len(sequences)},1-{end1 - 1}"
                    )
            # 不跨区域延伸
            else:
                # LSC 在前
                if (start2 - end1) <= (len(sequences) - end2 + start1):
                    result = (
                        f"LSC:{end2 + 1}-{len(sequences)}\tIRb:{start1 + 1}-{end1 - 1}\t"
                        f"SSC:{end1}-{start2 + 1}\tIRa:{start2 + 2}-{end2}"
                    )
                else:
                    result = (
                        f"LSC:{end1}-{start2 + 1}\tIRb:{start2 + 2}-{end2}\t"
                        f"SSC:{end2 + 1}-{len(sequences)}\tIRa:{start1 + 1}-{end1 - 1}"
                    )
        return result


def main():
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--input_file', help='fasta/GenBank format file', required=True)
    args = parser.parse_args()
    input_file = args.input_file
    try:
        result = check_files(input_file)
        if result:
            print(os.path.basename(input_file), end='\t')
            print(find_repeat_regions(result))
    except (ValueError, FileNotFoundError) as e:
        print(e)


if __name__ == '__main__':
    main()
