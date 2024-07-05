# CPStools Usage

## Install 

```sh
# python version >= 3.9
pip install biopython 
# 使用清华源安装 python version >= 3.9
pip install cpstools -i  https://pypi.tuna.tsinghua.edu.cn/simple
```

## Usage

```python
cpstools -h 
```

## subcommand

### gbcheck

gbcheck subcommand accepts two parameters:

 `-i, --test_file`

 `-r, --ref_file`

The parameter of -i/--test_file is required, and -r/--ref_file is optional.

If only -i/--TEST_FILE, the gbcheck.py will check the annotation of the input file. The script will check the gene name in each gene, CDS, tRNA and rRNA. For CDS, tRNA and rRNA, it will also check the 'product' label.

For CDS, the codon amino acids will be check, if start codon, stop codon and internal codon was wrong, it will print out the location of CDS.

```python
# self check
cpstools gbcheck -i input_file
# compara
cpstools gbcheck -i input_file -r ref_file
```



### info

info subcommand accepts one parameter:

 `-i, --input_file`

The parameter of -i/-input_file

This subcommand is used to statistic the gene types and gene numbers

```python
cpstools info -i input_file 
```



### IR

IR subcommand accepts one parameter:

 `-i, --input_file`

The parameter of -i/--input_file is required, and the input file can be genbank/fasta format. It will find the four regions in the chloroplast sequences, which the IR region will be identified if length surpass 1,000 bp.

```python
cpstools IR -i input_file
```



### Seq

The Seq subcommand can be used to adjust the start to the 1st bp in LSC region and the SSC forward. This analysis should combine with the co-linearity.

Seq subcommand accepts four parameter:

 `-d, --work_dir`

The parameter of -i/--work_dir is required, and the input directory is fasta format of all need adjusted sequences.

`-f, --info_file` 
The parameter offile of -f/--info_file is required, This file can be generated from IR subcommand.

-m {SSC,LSC,RP}, --mode {SSC,LSC,RP}

The parameter of -m is also required, And three choices was supported:
    `Mode: SSC for adjust_SSC_forward, 
          LSC foradjust_start_to_LSC, 
          RP for adjust sequence to reverse_complement`

```python
# adjust to LSC 1 bp start, the results are saved in LSC directory
cpstools Seq -d work_dir  -f info.txt -m LSC

# adjust SSC forward, the results are saved in SSC directory
cpstools Seq -d work_dir  -f info.txt -m SSC

# adjust to reversed complement and LSC 1 bp start, the results are saved in RP directory
cpstools Seq -d work_dir  -f info.txt -m RP
```



### Pi

The Pi subcommand accepts a required parameter `-d/--work_dir` and an optional parameter `-m/--mafft_path`. It automatically extracts shared gene regions and shared intergenic sequences from gb files in the work_dir, performs multiple sequence alignments, and calculates Pi values. Finally, it outputs the results sorted according to the chloroplast genome order.

```python
# If mafft in the environment variables
cpstools Pi -d work_dir

# if mafft not in the environment variables
cpstools Pi -d work_dir -m mafft_path
```



### RSCU

The RSCU subcommand accepts a required parameter `-d/--work_dir` and an optional parameter `-l/--filter_length`. It automatically extracts shared coding regions from GenBank files in the work_dir, removes duplicates, filters out short sequences, and calculates their RSCU values.

```python
# filter_length default is 300 bp
cpstools RSCU -d work_dir 
# filter_length is n bp
cpstools RSCU -d work_dir -l n
```

### SSRs

The SSRs subcommand accepts a required parameter `-i/--input_file` and an optional parameter `-k/--kmer_length`. It automatically identifies and locates SSRs in the GenBank files.

```python
# SSRs length, default is 10,6,5,4,4,4
cpstools SSRs -i input_file 
# SSRs length set to a specified length, 10,5,4,3,3,3
cpstools SSRs -i input_file -k 10,5,4,3,3,3
```

### converse

The converse subcommand accepts two required parameter `-d/--input_dir` and  `-m/--mode {fasta,mVISTA,tbl}`. It automatically converse genbank format file into fasta/mVISTA/tbl format.

```python
# converse genbank into fasta, the output in the fasta directory
cpstools converse -d input_dir -m fasta

# converse genbank into tbl, the output in the tbl directory
cpstools converse -d input_dir -m tbl

# converse genbank into mVISTA, the output in the mVISTA directory
cpstools converse -d input_dir -m mVISTA
```

### LSRs

The LSRs subcommand accepts one required parameter `-i/--input_file` . It automatically reads the output results from Reputer, sorts, and locates them.

```
cpstools LSRS -i input_file
```

### phy

The phy subcommand accepts two required parameter `-d/--input_dir` and  `-m/--mode {cds,pro}`. It automatically extracts shared CDS or protein sequences from GenBank files in the input_dir, sorts, and merges them.

```python
# to extract shared CDS sequences for phylogenetic analysis
cpstools phy -d input_dir -m cds

# to extract shared protein sequences for phylogenetic analysis
cpstools phy -d input_dir -m pro
```

