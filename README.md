# CPStools Usage

## Download:

```shell
git clone https://github.com/Xwb7533/CPStools.git
```

## Install dependencies

```sh
# python version >= 3.6
pip install biopython
```

## Usage

```python
python CPStools/bin/CPStools 
```

## subcommand

### gbcheck

gbcheck subcommand accepts two parameters:

 -i, --test_file

 -r, --ref_file

The parameter of -i/--test_file is required, and -r/--ref_file is optional.

If only -i/--TEST_FILE, the gbcheck.py will check the annotation of the input file. The script will check the gene name in each gene, CDS, tRNA and rRNA. For CDS, tRNA and rRNA, it will also check the 'product' label.

For CDS, the codon amino acids will be check, if start codon, stop codon and internal codon was wrong, it will print out the location of CDS.

```python
# self check
python CPStools/bin/CPStools gbcheck -i input_file
# compara
python CPStools/bin/CPStools gbcheck -i input_file -r ref_file
```

### info

info subcommand accepts two parameters:

 -i, --input_file

 -o, --output_file

The parameter of -i/-input_file and  -o/--output_file are required

This subcommand is used to statistic the gene types and gene numbers

```python
python CPStools/bin/CPStools info -i input_file -o aa.txt  
```

### IR

### Seq

### Pi_1

### Pi_2

### RSCU

### SSRs

### converse

### LSRs

### phy

