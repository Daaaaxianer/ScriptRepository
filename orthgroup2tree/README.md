# Orthgroups and phylogenetic tree work

All scripts support the `-h` parameter to view help.

## Dependency package

1. [Muscle 5](http://drive5.com/muscle/) not Muscle 3
2. [trimal V1.4](https://github.com/inab/trimal)
3. [iqtree2](https://github.com/iqtree/iqtree2)
4. Third-party python package:[Biopython](https://biopython.org/), pandas, numpy, subprocess, multiprocessing, argparse 

## Simple usage

**Function**: Process orthogroups and generate phylogenetic trees 

`orthogroup2tree.py [-h] [-i INDIR] [-o OUTDIR] [-p PROCESSES] [-e END] [-c COUNT] orthogroups_file`

```
positional arguments:

orthogroups_file      Path to the Orthogroups file

option arguments:

-h, --help                    show this help message and exit
-i INDIR, --indir INDIR       Directory containing genome files (default: genome)
-o OUTDIR, --outdir OUTDIR    Output directory name (default: out)
-p PROCESSES, --processes PROCESSES   Number of processes to use (default: 5)
-e END, --end END             Suffix for FASTA files (default: .pep)
-c COUNT, --count COUNT       Minimum count for row comparison (default: 5)

```

**example**
```
python orthogroup2tree.py orthogroups.test.txt
```