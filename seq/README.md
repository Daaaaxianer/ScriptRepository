# Seq Work Script
Sequence-related scripts are stored here.

All scripts support the `-h` parameter to view help.

## GffFormat.py
* Function 1: convert standard gff3 to custom gff (--gff3gff)
```
# -s is an optional parameter (default：species)
python GffFormat.py --gff3gff -g Vvinifera_145_Genoscope.12X.gene.gff3 -s vvi
```
* Function 2: read id file and extract bed for MG2C from custom gff (--id2bed)

The id file should have the same format as id.example.txt, with only one column: gene name
```
# -b is an optional parameter (default: out.bed)
python GffFormat.py --id2bed -i id.example.txt -f sl.4.0.all.gff -b sl.mg2c.bed
```

## SeqFormat.py
* Function 1: Translate nucleic acid coding sequence into amino acid sequence (--cdstrans)
```
# -p is an optional parameter (default: out.pep.fasta)
python SeqFormat.py --cdstrans -c sl.4.0.cds.fa -p out.pep.fasta
```
* Function 2: Calculate the gc content of nucleic acid sequences (--gc)
```
# -g is an optional parameter (default: out.gctable.txt)
python SeqFormat.py --gc -n sl.4.0.cds.fa -g out.gctable.txt
```
* Function 3: Extract the sequence of a given id (--extract)

The id file should have the same format as `id.example.txt`, with only one column: gene name
```
# -e is an optional parameter (default: out.extractedSeq.fasta)
python SeqFormat.py --extract -s sl.4.0.cds.fa -i id.example.txt -e out.extractedSeq.fasta
```
* Function 4: Cut the sequence according to the given location information (--cut)

The cut file should have the same format as `cut.example.txt`, which contains four columns: original id, start site, end site, new id.
```
# -m is an optional parameter (default: out.cutSeq.fasta)
python SeqFormat.py --cut -s sl.4.0.cds.fa -l cut.example.txt -m out.cutSeq.fasta
```



