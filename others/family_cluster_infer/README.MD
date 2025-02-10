# family cluster infer

## Based on the location of the family, determine the orthologous relationship of some clusters on the collinear block

using the -h parameter for help.

```
usage: family_cluster_infer.py [-h] --right_gff RIGHT_GFF --left_gff LEFT_GFF --right_id RIGHT_ID --left_id LEFT_ID --block_file BLOCK_FILE
                               [--window_size WINDOW_SIZE] [--min_cluster_size MIN_CLUSTER_SIZE] --output_dir OUTPUT_DIR

Process GFF and block files for clustering.

options:
  -h, --help            show this help message and exit
  --right_gff RIGHT_GFF
                        Path to Right.new.gff.txt file
  --left_gff LEFT_GFF   Path to Left.new.gff.txt file
  --right_id RIGHT_ID   Path to Right_id.txt file
  --left_id LEFT_ID     Path to Left_id.txt file
  --block_file BLOCK_FILE
                        Path to Left_Right.block.rr.txt file
  --window_size WINDOW_SIZE
                        Sliding window size in base pairs (default: 2e6)
  --min_cluster_size MIN_CLUSTER_SIZE
                        Minimum cluster size (default: 3)
  --output_dir OUTPUT_DIR
                        Directory path for output files

```

```
# example
python family_cluster_infer.py --left_gff Cume.new.gff --right_gff Bhhi.new.gff --left_id Cume_id.txt --right_id Bhhi_id.txt --block_file Cume_Bhhi.block.rr.txt --output_dir ./

```
