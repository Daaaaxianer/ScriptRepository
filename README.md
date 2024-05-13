# ScriptRepository

Sequence-related scripts are stored here.

All scripts support the -h parameter to view help.

## filterBlastBestPairs.py

Extraction of the best matching gene pairs from blast+ outfmt6

```
python filterBlastBestPairs.py -bl aaa.blast -o best.pairs.txt
```

## pyDotplot.py

Draw a genome dotplot map

```
# python pyDotplot.py query_name query_gff query_lens subject_name subject_gff subject_lens blastout pos
python pyDotplot.py Grape vv.12x.gff vv.12x.lens Tomato sl.4.0.all.gff sl.4.0.all.lens vv.sl.pep.1e-5.mts20.outfmt6 -p order
```

## ete3tree.py

Display multiple sequence alignments and tree files together
```
python ete3tree.py -h
```
