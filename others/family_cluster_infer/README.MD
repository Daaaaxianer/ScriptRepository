# family cluster infer

## Based on the location of the family, determine the orthologous relationship of some clusters on the collinear block

using the -h parameter for help.

```
usage: family_cluster_infer.py [-h] -odir ODIR [--window_size WINDOW_SIZE] [--min_cluster_size MIN_CLUSTER_SIZE]
                               [-lgff LGFF] [-rgff RGFF] [-lid LID] [-rid RID] [--family2cluster] [-lcluster LCLUSTER]
                               [-rcluster RCLUSTER] [-bkfile BKFILE] [--cluster2block]

Process GFF and block files for clustering.

options:
  -h, --help            show this help message and exit
  -odir ODIR            Directory path for output files
  --window_size WINDOW_SIZE
                        Sliding window size in base pairs (default: 2e6)
  --min_cluster_size MIN_CLUSTER_SIZE
                        Minimum cluster size (default: 3)

family2cluster:
  Arguments for family to cluster

  -lgff LGFF            Path to Left.new.gff.txt file
  -rgff RGFF            Path to Right.new.gff.txt file
  -lid LID              Path to Left_id.txt file
  -rid RID              Path to Right_id.txt file
  --family2cluster      Run family to cluster functionality

cluster2block:
  Arguments for cluster to block

  -lcluster LCLUSTER    Path to spLeftId2cluster.txt file
  -rcluster RCLUSTER    Path to spRightId2cluster.txt file
  -bkfile BKFILE        Path to Left_Right.block.rr.txt file
  --cluster2block       Run cluster to block functionality

```

## Function family2cluster：Calculate the family cluster based on physical location

```
# example
python family_cluster_infer.py --family2cluster -lgff Cume.new.gff -rgff Bhhi.new.gff -lid Cume_id.txt -rid Bhhi_id.txt -odir .
```

## Function cluster2block：Connecting clusters and collinear blocks

```
# example
python family_cluster_infer.py --cluster2block -lcluster spLeftId2cluster.txt -rcluster spRightId2cluster.txt -bkfile Cume_Bhhi.block.rr.txt -odir .
```



