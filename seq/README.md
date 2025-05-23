# Seq Work Script
Sequence-related scripts are stored here.

All scripts support the `-h` parameter to view help.

## GffFormat.py
* **Function 1**: convert standard gff3 to custom gff (--gff3gff)
```
# -s is an optional parameter (default：species)
python GffFormat.py --gff3gff -g Vvinifera_145_Genoscope.12X.gene.gff3 -s vvi
```
* **Function 2**: read id file and extract bed for MG2C from custom gff (--id2bed)

The id file should have the same format as `id.example.txt`, with only one column: gene name
```
# -b is an optional parameter (default: out.bed)
python GffFormat.py --id2bed -i id.example.txt -f sl.4.0.all.gff -b sl.mg2c.bed
```

* **Function 3**: Read gff filter file and extract id from custom gff (--filtergff) 

The gff filter file should be in the same format as `the T2Tcytobands.txt`, and the first column should be the same as the first column of the custom gff.

This function is the same as running `gff_filter.py` alone.
```
# -f is a reused parameter
python GffFormat.py --filtergff -f T2T.gff -ff T2Tcytobands.txt

# Equivalent to running gff_filter.py
python gff_filter.py -gff T2T.gff -fgff T2Tcytobands.txt
```

## SeqFormat.py
* **Function 1**: Translate nucleic acid coding sequence into amino acid sequence (--cdstrans)
```
# -p is an optional parameter (default: out.pep.fasta)
python SeqFormat.py --cdstrans -c sl.4.0.cds.fa -p out.pep.fasta
```
* **Function 2**: Calculate the gc content of nucleic acid sequences (--gc)
```
# -g is an optional parameter (default: out.gctable.txt)
python SeqFormat.py --gc -n sl.4.0.cds.fa -g out.gctable.txt
```
* **Function 3**: Extract the sequence of a given id (--extract)

The id file should have the same format as `id.example.txt`, with only one column: gene name
```
# -e is an optional parameter (default: out.extractedSeq.fasta)
python SeqFormat.py --extract -s sl.4.0.cds.fa -i id.example.txt -e out.extractedSeq.fasta
```
* **Function 4**: Cut the sequence according to the given location information (--cut)

The cut file should have the same format as `cut.example.txt`, which contains four columns: original id, start site, end site, new id.
```
# -m is an optional parameter (default: out.cutSeq.fasta)
python SeqFormat.py --cut -s sl.4.0.cds.fa -l cut.example.txt -m out.cutSeq.fasta
```

* **Function 5**: Split fatsa sequences by id (--splitfa)

This method can split the fasta sequence into strips of sequences according to the id name of the sequence.

```
# -odir is an optional parameter (default: split_out)
python SeqFormat.py --splitfa -fa sl.4.0.cds.fa
# Use-odir custom output directory
python SeqFormat.py --splitfa -fa sl.4.0.cds.fa -odir out

# Equivalent to running split_fasta.py
python split_fasta.py -fa sl.4.0.cds.fa -odir iout

```

## del_seq.py

```
usage: del_seq.py [-h] [--min_length MIN_LENGTH] [--exclude_genes EXCLUDE_GENES] [--cds_output CDS_OUTPUT]
                  [--pep_output PEP_OUTPUT]
                  cds_file pep_file

过滤CDS和PEP序列，可按长度过滤或排除指定基因

positional arguments:
  cds_file              输入的CDS序列文件(FASTA格式)
  pep_file              输入的PEP序列文件(FASTA格式)

options:
  -h, --help            show this help message and exit
  --min_length MIN_LENGTH
                        最小长度阈值(默认200)，当不使用--exclude_genes时有效
  --exclude_genes EXCLUDE_GENES
                        包含要排除基因列表的文件(单列基因名)，如果使用此参数，将忽略--min_length
  --cds_output CDS_OUTPUT
                        输出的CDS序列文件名(默认filtered_cds.fasta)
  --pep_output PEP_OUTPUT
                        输出的PEP序列文件名(默认filtered_pep.fasta)
```

* **Function 1**: delete sequence sequences by length (--min_length)

```
## Delete sequences with a length of less than 200
python del_seq.py cds_file pep_file --min_length 200

## Delete sequences with a length of less than 100 and specify the output filename
python del_seq.py cds_file pep_file --min_length 100 --cds_output xxx.less100.cds.fasta --pep_output xxx.less100.pep.fasta
```

* **Function 2**: delete sequence sequences by id (--exclude_genes)

If this parameter is used, the --min_length parameter will be ignored.
```
## Remove the genes listed with their IDs from the file.
python del_seq.py cds_file pep_file --exclude_genes xxx.id.txt

## Remove the genes listed with their IDs from the file and specify the output filename.
python del_seq.py cds_file pep_file --exclude_genes xxx.id.txt --cds_output xxx.removeid.cds.fasta --pep_output xxx.removeid.pep.fasta
```
  






