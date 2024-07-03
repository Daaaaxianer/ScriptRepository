# Orthgroups and phylogenetic tree work

All scripts support the `-h` parameter to view help.

## Dependency package

1. (Muscle 5)[http://drive5.com/muscle/] not Muscle 3
2. (trimal V1.4)[https://github.com/inab/trimal]
3. (iqtree2)[https://github.com/iqtree/iqtree2]
4. python package:(Biopython)[https://biopython.org/], pandas, numpy, subprocess, multiprocessing, argparse 

## Simple usage

* **Function **: Process orthogroups and generate phylogenetic trees (--gff3gff)
```
# -s is an optional parameter (default：species)
python GffFormat.py --gff3gff -g Vvinifera_145_Genoscope.12X.gene.gff3 -s vvi
```