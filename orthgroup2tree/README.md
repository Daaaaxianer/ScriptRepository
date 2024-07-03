# Orthgroups and phylogenetic tree work

All scripts support the `-h` parameter to view help.

## Dependency package

1. [Muscle 5](http://drive5.com/muscle/) not Muscle 3
2. [trimal V1.4](https://github.com/inab/trimal)
3. [iqtree2](https://github.com/iqtree/iqtree2)
4. Third-party python package:[Biopython](https://biopython.org/), pandas, numpy, subprocess, multiprocessing, argparse 

## Simple usage

**Function **: Process orthogroups and generate phylogenetic trees 

`orthogroup2tree.py [-h] [-i INDIR] [-o OUTDIR] [-p PROCESSES] [-e END] [-c COUNT] orthogroups_file`

positional arguments:

orthogroups_file      Path to the Orthogroups file

```
# -i Directory containing genome files (default: genome)
# -o Output directory name (default: out)
# -p Number of processes to use (default: 5)
# -e Suffix for FASTA files (default: .pep)
# -c Minimum count for row comparison (default: 5)

python GffFormat.py --gff3gff -g Vvinifera_145_Genoscope.12X.gene.gff3 -s vvi
```