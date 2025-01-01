# codonBias

Using [CodonU](https://github.com/SouradiptoC/CodonU) to Analyze the Preference of Codons

# Dependency

We need a lot of Python packages to perform the following operations

```
pip install -r codon_requirements.txt
```

All scripts support the `-h` parameter to view help.

### Step One: Generate counts and reports 

**Usage**
```
python 1.codonBias_report.py <nucleotide_file> <protein_file> <output_directory>
```
**Example**
```
python 1.codonBias_report.py 2OGD.in.cc.cds.fasta 2OGD.in.cc.protein.fasta 2OGD.in.cc
```
### Step Two: Generate plots 

**Usage**
```
python 2.codonBias_report.py <nucleotide_file> <protein_file> <output_directory>
```
