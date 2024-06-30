import argparse
from Bio import SeqIO
from collections import defaultdict


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Compare gene counts and difference in two GenBank files.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com
        """

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


def check_gene(feature):
    if 'gene' not in feature.qualifiers:
        print(f"The gene is missing a name at location {feature.location}\n{'-' * 40}")


def check_CDS(feature, extract_seq):
    if 'gene' in feature.qualifiers:
        gene_name = feature.qualifiers['gene'][0]
        if len(extract_seq) % 3 != 0:
            print(f"The CDS: {gene_name} at {feature.location} length {len(extract_seq)}, "
                  f"is not a multiple of 3!\n{'-' * 40}")
        elif str(extract_seq)[:3] != "ATG":
            check_non_atg_CDS(feature, extract_seq, gene_name)
    else:
        print(f"The CDS is missing a gene name at location {feature.location}\n{'-' * 40}")


def check_tRNA_rRNA(feature):
    if 'gene' not in feature.qualifiers:
        print(f"The {feature.type} is missing a gene name at location {feature.location}\n{'-' * 40}")
    if 'product' not in feature.qualifiers:
        print(f"The {feature.type} is missing a product at location {feature.location}\n{'-' * 40}")


def check_non_atg_CDS(feature, extract_seq, gene_name):
    my_pro = extract_seq.translate(table=11)
    if "exception" in feature.qualifiers:
        except_type = feature.qualifiers["exception"][0]
        print(f"The CDS: {gene_name} at {feature.location} does not start with 'ATG' "
              f"and is annotated as {except_type}\n{'-' * 40}")
        if not my_pro[-1] == "*":
            print(f"The CDS: {gene_name} at {feature.location} stop codon is incorrect!\n{'-' * 40}")
        elif "*" in my_pro[:-1]:
            print(f"The CDS: {gene_name} at {feature.location} has an internal stop codon!\n{'-' * 40}")
    else:
        print(f"The CDS: {gene_name} at {feature.location} does not start with 'ATG', please verify!\n{'-' * 40}")
        if not my_pro[-1] == "*":
            print(f"The CDS: {gene_name} at {feature.location} stop codon is incorrect!\n{'-' * 40}")
        elif "*" in my_pro[:-1]:
            print(f"The CDS: {gene_name} at {feature.location} has an internal stop codon!\n{'-' * 40}")


def check(input_file):
    records = parse_genbank_file(input_file)
    if records:
        for rec in records:
            for feature in rec.features:
                if feature.type == 'gene':
                    check_gene(feature)
                elif feature.type == 'CDS':
                    check_CDS(feature, feature.extract(rec.seq))
                elif feature.type in ['rRNA', 'tRNA']:
                    check_tRNA_rRNA(feature)
    else:
        print(f"{input_file} is not a valid Genbank format file")


def compare_genbank_files(ref_file, test_file):
    ref_dict, test_dict = defaultdict(int), defaultdict(int)

    # ref file
    ref_records = parse_genbank_file(ref_file)
    if ref_records is None:
        return

    for rec in ref_records:
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'rRNA' or feature.type == 'tRNA':
                if feature.qualifiers.get('gene'):
                    gene_name = feature.qualifiers.get('gene')[0]
                    ref_dict[gene_name] += 1
                else:
                    gene_loc = feature.location
                    gene_start = int(gene_loc.start) + 1
                    gene_end = int(gene_loc.end)
                    print(f"[{gene_start, gene_end}]  in reference gb file has no gene name!")

    # test file
    test_records = parse_genbank_file(test_file)
    if test_records is None:
        return
    print('\n' + '-'*80)
    for rec in test_records:
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'rRNA' or feature.type == 'tRNA':
                if feature.qualifiers.get('gene'):
                    gene_name = feature.qualifiers.get('gene')[0]
                    test_dict[gene_name] += 1
                else:
                    gene_loc = feature.location
                    gene_start = int(gene_loc.start) + 1
                    gene_end = int(gene_loc.end)
                    print(f"[{gene_start, gene_end}]  in test gb file has no gene name!")

    print(f"{ref_file} file has {len(ref_dict)} genes!")
    print(f"{test_file} has {len(test_dict)} genes!")

    # Comparative
    for gene in ref_dict.keys():
        if ref_dict[gene] != test_dict[gene]:
            print(f"{gene} in {ref_file} number is {ref_dict[gene]} \
                but in {test_file} is {test_dict[gene]}")
    result_dict = {key: value for key, value in test_dict.items() if key not in ref_dict}
    print(f"Un annotated Genes in {ref_file} is {result_dict}")
    print('-'*80 + '\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-r', '--ref_file', help='reference GenBank file')
    parser.add_argument('-i', '--test_file', help='testing GenBank file', required=True)
    args = parser.parse_args()
    if args.test_file and not args.ref_file:
        print(args.test_file)
        check(args.test_file)
    elif args.test_file and args.ref_file:
        compare_genbank_files(args.ref_file, args.test_file)
    else:
        parser.print_help()
