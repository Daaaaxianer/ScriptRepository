# Others

## TE_classification_and_counting

Categorize and count the results of TE.

```
# default
python te_classification_and_counting.py GWHBCHF00000000.genome.fasta.out

# "Optional parameter -o output_path"
python te_classification_and_counting.py GWHBCHF00000000.genome.fasta.out -o ./
```

## interproscan_classification_and_counting

Categorize and count the results of interproscan.

```
python interproscan_classification_and_counting.py interproscan.result.tsv interproscan.result.count.txt
```

## WGDI file to JCVI file

Generate input files for JCVI using the input or output results of WGDI



