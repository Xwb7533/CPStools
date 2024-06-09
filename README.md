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

IR subcommand accepts one parameter:

 -i, --input_file

The parameter of -i/--input_file is required, and the input file can be genbank/fasta format. It will find the four regions in the chloroplast sequences, which the IR region will be identified if length surpass 1,000 bp.

```python
python CPStools/bin/CPStools IR -i input_file
```

### Seq

The Seq subcommand can be used to adjust the start to the 1st bp in LSC region and the SSC forward. This analysis should combine with the co-linearity.

Seq subcommand accepts four parameter:

 -i, --work_dir

The parameter of -i/--work_dir is required, and the input directory is fasta format of all need adjusted sequences.

 -o, --save_dir

The parameter of -o/--save_dir is required, The directry of adjusted sequences.

-f, --info_file 
The parameter offile of -f/--info_file is required, This file can be generated from IR subcommand.

-m {SSC,LSC,RP}, --mode {SSC,LSC,RP}

The parameter of -m is also required, And three choices was supported:
    Mode: SSC for adjust_SSC_forward, 
          LSC foradjust_start_to_LSC, 
          RP for adjust sequence to reverse_complement

```python
python CPStools/bin/CPStools Seq -i input_dir -o output_dir -f info,txt -m LSC/SSC/RP
```

### Pi_1

### Pi_2

### RSCU

### SSRs

### converse

### LSRs

### phy

