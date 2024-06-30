from Bio import SeqIO
import os, sys
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        To get information and intron numbers from genbank files.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""

        return f"{description}\n\n{help_text}"

# information table

def information_table(input_file, output_file):
    subunits_of_photosystem_I = {}
    Subunits_of_photosystem_II = {}
    Subunits_of_NADH_dehydrogenase = {}
    Subunits_of_cytochrome_b_f_complex = {}
    Large_subunit_of_rubisco = {}
    Subunits_reductase = {}
    Subunits_of_RNA_polymerase = {}
    Proteins_of_large_ribosomal_subunit = {}
    Proteins_of_small_ribosomal_subunit = {}
    Subunits_of_ATP_synthase = {}
    Maturase = {}
    Protease = {}
    others = {}
    ycf_gene = {}
    ORFS = {}
    trnA_gene = {}
    rrnA_gene = {}
    Envelope_membrane_protein = {}
    Acetyl_CoA, Translation = {}, {}
    cytochrome_synthesis = {}

    for rec in SeqIO.parse(input_file, 'genbank'):
        for feature in rec.features:
            if feature.type == 'gene':
                gene_name = feature.qualifiers['gene'][0]
                if gene_name.lower().startswith('rps'):
                    if gene_name not in Proteins_of_small_ribosomal_subunit:
                        Proteins_of_small_ribosomal_subunit[gene_name] = 1
                    else:
                        Proteins_of_small_ribosomal_subunit[gene_name] += 1
                elif gene_name.lower().startswith('psa'):
                    if gene_name not in subunits_of_photosystem_I.keys():
                        subunits_of_photosystem_I[gene_name] = 1
                    else:
                        subunits_of_photosystem_I[gene_name] += 1
                elif gene_name.lower().startswith('psb'):
                    if gene_name not in Subunits_of_photosystem_II.keys():
                        Subunits_of_photosystem_II[gene_name] = 1
                    else:
                        Subunits_of_photosystem_II[gene_name] += 1
                elif gene_name.lower().startswith('ndh'):
                    if gene_name not in Subunits_of_NADH_dehydrogenase.keys():
                        Subunits_of_NADH_dehydrogenase[gene_name] = 1
                    else:
                        Subunits_of_NADH_dehydrogenase[gene_name] += 1
                elif gene_name.lower().startswith('pet'):
                    if gene_name not in Subunits_of_cytochrome_b_f_complex.keys():
                        Subunits_of_cytochrome_b_f_complex[gene_name] = 1
                    else:
                        Subunits_of_cytochrome_b_f_complex[gene_name] += 1
                elif gene_name.lower().startswith('rbcl'):
                    if gene_name not in Large_subunit_of_rubisco.keys():
                        Large_subunit_of_rubisco[gene_name] = 1
                    else:
                        Large_subunit_of_rubisco[gene_name] += 1
                elif gene_name.lower().startswith('rpo'):
                    if gene_name not in Subunits_of_RNA_polymerase.keys():
                        Subunits_of_RNA_polymerase[gene_name] = 1
                    else:
                        Subunits_of_RNA_polymerase[gene_name] += 1
                elif gene_name.lower().startswith('rpl'):
                    if gene_name not in Proteins_of_large_ribosomal_subunit.keys():
                        Proteins_of_large_ribosomal_subunit[gene_name] = 1
                    else:
                        Proteins_of_large_ribosomal_subunit[gene_name] += 1
                elif gene_name.lower().startswith('atp'):
                    if gene_name not in Subunits_of_ATP_synthase.keys():
                        Subunits_of_ATP_synthase[gene_name] = 1
                    else:
                        Subunits_of_ATP_synthase[gene_name] += 1
                elif gene_name.lower().startswith('matk'):
                    if gene_name not in Maturase.keys():
                        Maturase[gene_name] = 1
                    else:
                        Maturase[gene_name] += 1
                elif gene_name.lower().startswith('clpp'):
                    if gene_name not in Protease.keys():
                        Protease[gene_name] = 1
                    else:
                        Protease[gene_name] += 1
                elif gene_name.lower().startswith('ycf'):
                    if gene_name not in ycf_gene.keys():
                        ycf_gene[gene_name] = 1
                    else:
                        ycf_gene[gene_name] += 1
                elif gene_name.lower().startswith('orf'):
                    if gene_name not in ORFS.keys():
                        ORFS[gene_name] = 1
                    else:
                        ORFS[gene_name] += 1
                elif gene_name.lower().startswith('trn'):
                    if gene_name not in trnA_gene.keys():
                        trnA_gene[gene_name] = 1
                    else:
                        trnA_gene[gene_name] += 1
                elif gene_name.lower().startswith('rrn'):
                    if gene_name not in rrnA_gene.keys():
                        rrnA_gene[gene_name] = 1
                    else:
                        rrnA_gene[gene_name] += 1
                elif gene_name.lower().startswith('ch'):
                    if gene_name not in Subunits_reductase.keys():
                        Subunits_reductase[gene_name] = 1
                    else:
                        Subunits_reductase[gene_name] += 1
                elif gene_name.lower().startswith('ce'):
                    if gene_name not in Envelope_membrane_protein.keys():
                        Envelope_membrane_protein[gene_name] = 1
                    else:
                        Envelope_membrane_protein[gene_name] += 1
                elif gene_name.lower().startswith('acc'):
                    if gene_name not in Acetyl_CoA.keys():
                        Acetyl_CoA[gene_name] = 1
                    else:
                        Acetyl_CoA[gene_name] += 1
                elif gene_name.lower().startswith('ccs'):
                    if gene_name not in cytochrome_synthesis.keys():
                        cytochrome_synthesis[gene_name] = 1
                    else:
                        cytochrome_synthesis[gene_name] += 1
                elif gene_name.lower().startswith('inf'):
                    if gene_name not in Translation.keys():
                        Translation[gene_name] = 1
                    else:
                        Translation[gene_name] += 1
                else:
                    if gene_name not in others.keys():
                        others[gene_name] = 1
                    else:
                        others[gene_name] += 1
    with open(output_file, 'w') as ff:
        gene_1, gene_2, gene_3, gene_4, gene_5, gene_6, gene_7, gene_8 = [], [], [], [], [], [], [], []
        gene_9, gene_10, gene_11, gene_12, gene_13, gene_14, gene_15, gene_16 = [], [], [], [], [], [], [], []
        gene_17, gene_18, gene_19, gene_20 = [], [], [], []
        ff.write('Category\tGene group\tGene name\n')
        for key, value in subunits_of_photosystem_I.items():
            if value == 1:
                gene_1.append(key)
            else:
                gene_1.append(f"{key}({value})")
            gene_1_ = ','.join(i for i in gene_1)
        ff.write(f'Photosynthesis\tSubunits of photosystem I\t{gene_1_}\n')
        for key, value in Subunits_of_photosystem_II.items():
            if value == 1:
                gene_2.append(key)
            else:
                gene_2.append(f"{key}({value})")
            gene_2_ = ','.join(i for i in gene_2)
        ff.write(f'\tSubunits of photosystem II\t{gene_2_}\n')
        for key, value in Subunits_of_NADH_dehydrogenase.items():
            if value == 1:
                gene_3.append(key)
            else:
                gene_3.append(f"{key}({value})")
            gene_3_ = ','.join(i for i in gene_3)
        ff.write(f'\tSubunits of NADH dehydrogenase\t{gene_3_}\n')
        for key, value in Subunits_of_cytochrome_b_f_complex.items():
            if value == 1:
                gene_4.append(key)
            else:
                gene_4.append(f"{key}({value})")
            gene_4_ = ','.join(i for i in gene_4)
        ff.write(f'\tSubunits of cytochrome b/f complex\t{gene_4_}\n')
        for key, value in Subunits_of_ATP_synthase.items():
            if value == 1:
                gene_6.append(key)
            else:
                gene_6.append(f"{key}({value})")
            gene_6_ = ','.join(i for i in gene_6)
        ff.write(f'\tSubunits of ATP synthase\t{gene_6_}\n')
        for key, value in Large_subunit_of_rubisco.items():
            if value == 1:
                gene_5.append(key)
            else:
                gene_5.append(f"{key}({value})")
            gene_5_ = ','.join(i for i in gene_5)
        ff.write(f'\tLarge subunit of rubisco\t{gene_5_}\n')
        if Subunits_reductase:
            for key, value in Subunits_reductase.items():
                if value == 1:
                    gene_7.append(key)
                else:
                    gene_7.append(f"{key}({value})")
                gene_7_ = ','.join(i for i in gene_7)
            ff.write(f'\tSubunits photochlorophyllide reductase\t{gene_7_}\n')
        for key, value in Proteins_of_large_ribosomal_subunit.items():
            if value == 1:
                gene_8.append(key)
            else:
                gene_8.append(f"{key}({value})")
            gene_8_ = ','.join(i for i in gene_8)
        ff.write(f'Self-replication\tProteins of large ribosomal subunit\t{gene_8_}\n')
        for key, value in Proteins_of_small_ribosomal_subunit.items():
            if value == 1:
                gene_9.append(key)
            else:
                gene_9.append(f"{key}({value})")
            gene_9_ = ','.join(i for i in gene_9)
        ff.write(f'\tProteins of small ribosomal subunit\t{gene_9_}\n')
        for key, value in Subunits_of_RNA_polymerase.items():
            if value == 1:
                gene_10.append(key)
            else:
                gene_10.append(f"{key}({value})")
            gene_10_ = ','.join(i for i in gene_10)
        ff.write(f'\tSubunits of RNA polymerase\t{gene_10_}\n')
        for key, value in rrnA_gene.items():
            if value == 1:
                gene_11.append(key)
            else:
                gene_11.append(f"{key}({value})")
            gene_11_ = ','.join(i for i in gene_11)
        ff.write(f'\tRibosomal RNAs\t{gene_11_}\n')
        for key, value in trnA_gene.items():
            if value == 1:
                gene_12.append(key)
            else:
                gene_12.append(f"{key}({value})")
            gene_12_ = ','.join(i for i in gene_12)
        ff.write(f'\tTransfer RNAs\t{gene_12_}\n')
        for key, value in Maturase.items():
            if value == 1:
                gene_13.append(key)
            else:
                gene_13.append(f"{key}({value})")
            gene_13_ = ','.join(i for i in gene_13)
        ff.write(f'Other genes\tMaturase\t{gene_13_}\n')
        if Protease:
            for key, value in Protease.items():
                if value == 1:
                    gene_14.append(key)
                else:
                    gene_14.append(f"{key}({value})")
                gene_14_ = ','.join(i for i in gene_14)
            ff.write(f'\tProtease\t{gene_14_}\n')
        if Envelope_membrane_protein:
            for key, value in Envelope_membrane_protein.items():
                if value == 1:
                    gene_15.append(key)
                else:
                    gene_15.append(f"{key}({value})")
                gene_15_ = ','.join(i for i in gene_15)
            ff.write(f'\tEnvelope membrane protein\t{gene_15_}\n')
        if Acetyl_CoA:
            for key, value in Acetyl_CoA.items():
                if value == 1:
                    gene_16.append(key)
                else:
                    gene_16.append(f"{key}({value})")
                gene_16_ = ','.join(i for i in gene_16)
            ff.write(f'\tAcetyl-CoA carboxylase\t{gene_16_}\n')
        if cytochrome_synthesis:
            for key, value in cytochrome_synthesis.items():
                if value == 1:
                    gene_17.append(key)
                else:
                    gene_17.append(f"{key}({value})")
                gene_17_ = ','.join(i for i in gene_17)
            ff.write(f'\tc-type cytochrome synthesis gene\t{gene_17_}\n')
        if Translation:
            for key, value in Translation.items():
                if value == 1:
                    gene_19.append(key)
                else:
                    gene_19.append(f"{key}({value})")
                gene_19_ = ','.join(i for i in gene_19)
            ff.write(f'\tTranslation initiation factor\t{gene_19_}\n')
        if others:
            for key, value in others.items():
                if value == 1:
                    gene_18.append(key)
                else:
                    gene_18.append(f"{key}({value})")
                gene_18_ = ','.join(i for i in gene_18)
            ff.write(f'\tother\t{gene_18_}\n')


def intron_find(input_file):
    two_exon, three_exon, other_exon = [], [], []
    for record in SeqIO.parse(input_file, 'genbank'):
        for feat in record.features:
            if feat.type == 'CDS' or feat.type == 'tRNA' or feat.type == 'rRNA':
                if len(feat.location.parts) == 2:
                    if feat.qualifiers['gene'][0] not in two_exon:
                        two_exon.append(feat.qualifiers['gene'][0])
                elif len(feat.location.parts) == 3:
                    if feat.qualifiers['gene'][0] not in three_exon:
                        three_exon.append(feat.qualifiers['gene'][0])
                elif len(feat.location.parts) > 3:
                    if feat.qualifiers['gene'][0] not in three_exon:
                        other_exon.append(feat.qualifiers['gene'][0])
    print(f"Two exons {two_exon} \n")
    print(f"Three exons {three_exon}\n")
    if other_exon:
        print(f"More than three exons {other_exon}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--input_file", help="Input genbank format file")
    # parser.add_argument("-o", "--output_file", help="Output TSV file")
    args = parser.parse_args()
    if args.input_file:
        if args.input_file.endswith('gb') or args.input_file.endswith('gbk'):
            input_file_path = os.path.abspath(args.input_file)
            output_file_path = input_file_path.split('.')[0] + '.tsv'
            information_table(args.input_file, output_file_path)
            intron_find(args.input_file)
            # abs_path = os.path.abspath(args.output_file)
            print(f"The statistic results has been saved into {output_file_path}.")
        else:
            print(f"please input genbank format files and endwith 'gb' or 'gbk'! ")
    else:
        parser.print_help()




