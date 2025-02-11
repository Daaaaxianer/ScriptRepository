# family cluster infer

## Based on the location of the family, determine the orthologous relationship of some clusters on the collinear block

using the -h parameter for help.

```
usage: family_cluster_infer.py [-h] [-odir ODIR] [-gff GFF] [-idfile IDFILE] [-window_size WINDOW_SIZE]
                               [-min_cluster_size MIN_CLUSTER_SIZE] [-sp SP] [--family2cluster] [-lcluster LCLUSTER]
                               [-rcluster RCLUSTER] [-bkfile BKFILE] [--cluster2block]

Process GFF and block files for clustering.

options:
  -h, --help            show this help message and exit
  -odir ODIR            Directory path for output files (default: current directory)

family2cluster:
  Arguments for family to cluster

  -gff GFF              Path to the GFF file (e.g., Left.new.gff.txt)
  -idfile IDFILE        Path to the ID file (e.g., Left_id.txt)
  -window_size WINDOW_SIZE
                        Sliding window size in base pairs (default: 2e6)
  -min_cluster_size MIN_CLUSTER_SIZE
                        Minimum cluster size (default: 3)
  -sp SP                Output file prefix (default: sp)
  --family2cluster      Run family to cluster functionality

cluster2block:
  Arguments for cluster to block

  -lcluster LCLUSTER    Path to spLeftId2cluster.txt file
  -rcluster RCLUSTER    Path to spRightId2cluster.txt file
  -bkfile BKFILE        Path to Left_Right.block.rr.txt file
  --cluster2block       Run cluster to block functionality

```

## Function family2cluster：Calculate the family cluster based on physical location

**Input**：`Left.gff`, `Left.id.txt`

**Ouput**: `sp.id2cluster.txt`

```
# default example1
python family_cluster_infer.py --family2cluster -gff Cume.new.gff -idfile Cume_id.txt -sp cume

# default example2
python family_cluster_infer.py --family2cluster -gff Bhhi.new.gff -idfile Bhhi_id.txt -sp bhhi

# set -window_size,-min_cluster_size and -odir
python family_cluster_infer.py --family2cluster -gff Cume.new.gff -idfile Cume_id.txt -sp cume -window_size 3e6 -min_cluster_size 4 -odir ./output
```

## Function cluster2block：Connecting clusters and collinear blocks

**Input**：`spLeft.Id2cluster.txt`, `spRight.Id2cluster.txt` and `spLeft_spRight.block.rr.txt`

**Output**:`Left_Right.block.rr.txt`

```
# default example
python family_cluster_infer.py --cluster2block -lcluster cume.id2cluster.txt -rcluster bhhi.id2cluster.txt -bkfile Cume_Bhhi.block.rr.txt

# set -odir
python family_cluster_infer.py --cluster2block -lcluster cume.id2cluster.txt -rcluster bhhi.id2cluster.txt -bkfile Cume_Bhhi.block.rr.txt -odir ./output
```



