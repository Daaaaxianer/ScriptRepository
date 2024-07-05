# Mutiple alignment and Evolution tree
## Dependency package
* [Muscle 5](https://www.drive5.com/muscle5/)
* [trimAl V1.4]((https://github.com/inab/trimal)/

  It is recommended to compile with the latest version of the code on github
* [iqtree2](https://github.com/iqtree/iqtree2)
## muscle2iqtree.py
```
# muscle5+iqtree2 (linux and windows)
python muscle2iqtree.py xxx.protein.fasta
```
## muscle2trimal2iqtree.py (faster and recommended)
```
# muscle5+trimal1.4+iqtree2 (linux)
python muscle2trimal2iqtree.py xxx.cds.fasta
```
