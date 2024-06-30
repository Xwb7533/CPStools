from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Trans Geneious Annotated file into GB format.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com
        """

        return f"{description}\n\n{help_text}"

sort_order = {'gene': 0, 'CDS': 1, 'tRNA': 2, 'rRNA': 3}


def filter_qualifiers(feature, rec):
    # List of qualifiers to keep

    keys_to_keep = ['gene', 'location', 'exception', 'trans_splicing', 'codon_start',
                    'transl_table', 'product', 'translation', 'organism', 'organelle', 'mol_type']
    # Create a new feature with only the desired qualifiers
    filtered_qualifiers = {key: feature.qualifiers[key] for key in keys_to_keep if key in feature.qualifiers}
    if feature.type == 'CDS' and 'translation' not in filtered_qualifiers.keys():
        filtered_qualifiers['translation'] = feature.extract(rec.seq).translate(table=11)[:-1]
    return SeqFeature(location=feature.location, type=feature.type, qualifiers=filtered_qualifiers)


def check(input_file, output_file):
    with open(output_file, 'w') as output_handle:
        for rec in SeqIO.parse(input_file, 'gb'):
            # Skip the 'source' feature
            sorted_features = sorted((feature for feature in rec.features if feature.type != 'source'),
                                     key=lambda feature: (int(feature.location.start), int(feature.location.end), # Ensure position is integer
                                                          sort_order.get(feature.type, 4)))

            # Include the 'source' feature at the beginning
            source_feature = next((feature for feature in rec.features if feature.type == 'source'), None)
            if source_feature:
                sorted_features.insert(0, source_feature)

            # Filter qualifiers for each feature
            sorted_features_with_filtered_qualifiers = [filter_qualifiers(feature, rec) for feature in sorted_features]

            # Create a new SeqRecord with the sorted features
            new_rec = SeqRecord(rec.seq, id=rec.id, name=rec.name, description=rec.description,
                                annotations=rec.annotations, features=sorted_features_with_filtered_qualifiers)

            # Write the modified SeqRecord to the output file
            SeqIO.write(new_rec, output_handle, 'genbank')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--input_file', help='Geneious format file', required=True)
    parser.add_argument('-o', '--output_file', help='GenBank format file', required=True)
    args = parser.parse_args()
    if args.input_file and args.output_file:
        check(args.input_file, args.output_file)
    else:
        parser.print_help()
