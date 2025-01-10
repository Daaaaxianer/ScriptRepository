# Others

## Function One：TE_classification_and_counting

Categorize and count the results of TE.

```
# default
python te_classification_and_counting.py GWHBCHF00000000.genome.fasta.out

# "Optional parameter -o output_path"
python te_classification_and_counting.py GWHBCHF00000000.genome.fasta.out -o ./
```

## Function Two：interproscan_classification_and_counting

Categorize and count the results of interproscan.

```
python interproscan_classification_and_counting.py interproscan.result.tsv interproscan.result.count.txt
```

## Function Three：WGDI file to JCVI file

Generate input files for JCVI using the input or output results of WGDI

```
python wgdi2jcvi.py -h
```

### Usage 1 :wgdi gff to jcvi bed
``` 
python wgdi2jcvi.py gff2bed [-h] -gff input_gff_file -bed output_bed_file

# -gff input_gff_file 
# -bed output_bed_file
```
```
# example
python wgdi2jcvi.py gff2bed -gff coca_gff -bed coca.bed
```

### Usage 2 :blast to jcvi anchor




