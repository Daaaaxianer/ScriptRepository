# Seq Work Script
Sequence-related scripts are stored here.

All scripts support the `-h` parameter to view help.

## GffFormat.py
* Function 1: convert standard gff3 to custom gff
```
# example
# -s is an optional parameter (default：species)
python GffFormat.py --gff3gff -g Vvinifera_145_Genoscope.12X.gene.gff3 -s vvi
```
* Function 2: read id file and extract bed for MG2C from custom gff
```
# example
# -b is an optional parameter (default: out.bed)
python GffFormat.py --id2bed -i id.example.txt -f sl.4.0.all.gff -b sl.mg2c.bed
```
## SeqFormat.py
* Function 1: convert standard gff3 to custom gff
