# -*- coding:UTF-8 -*-
# FileName  :family_cluster_infer.py
# Time      :2025/2/10 
# Author    :xian

import pandas as pd
import re
import argparse
import os  # 引入 os 模块


# 添加命令行参数解析
def parse_args():
    parser = argparse.ArgumentParser(description="Process GFF and block files for clustering.")
    
    # 公共参数
    parser.add_argument("-odir", type=str, required=True, help="Directory path for output files")
    parser.add_argument("--window_size", type=float, default=2e6, help="Sliding window size in base pairs (default: 2e6)")
    parser.add_argument("--min_cluster_size", type=int, default=3, help="Minimum cluster size (default: 3)")

    # 定义 --family2cluster 功能所需的参数
    family2cluster = parser.add_argument_group('family2cluster', 'Arguments for family to cluster')
    family2cluster.add_argument("-lgff", type=str, required=False, help="Path to Left.new.gff.txt file")
    family2cluster.add_argument("-rgff", type=str, required=False, help="Path to Right.new.gff.txt file")
    family2cluster.add_argument("-lid", type=str, required=False, help="Path to Left_id.txt file")
    family2cluster.add_argument("-rid", type=str, required=False, help="Path to Right_id.txt file")
    family2cluster.add_argument("--family2cluster", action="store_true", help="Run family to cluster functionality")

    # 定义 --cluster2block 功能所需的参数
    cluster2block = parser.add_argument_group('cluster2block', 'Arguments for cluster to block')
    cluster2block.add_argument("-lcluster", type=str, required=False, help="Path to spLeftId2cluster.txt file")
    cluster2block.add_argument("-rcluster", type=str, required=False, help="Path to spRightId2cluster.txt file")
    cluster2block.add_argument("-bkfile", type=str, required=False, help="Path to Left_Right.block.rr.txt file")
    cluster2block.add_argument("--cluster2block", action="store_true", help="Run cluster to block functionality")

    return parser.parse_args()


# 主函数
def main():
    args = parse_args()

    # 检查并创建输出目录
    if not os.path.exists(args.odir):
        os.makedirs(args.odir)
        print(f"Output directory '{args.odir}' created.")

    if args.family2cluster:
        print("Running family2cluster functionality...")
        run_family2cluster(args)
    elif args.cluster2block:
        print("Running cluster2block functionality...")
        run_cluster2block(args)
    else:
        print("Please specify either --family2cluster or --cluster2block.")


# 运行 family2cluster 功能
def run_family2cluster(args):
    if not all([args.lgff, args.rgff, args.lid, args.rid]):
        print("Error: For --family2cluster, you must provide -lgff, -rgff, -lid, and -rid.")
        return

    # 读取 GFF 文件
    right_gff_cols = ['chr', 'start', 'end', 'strand', 'oldid', 'id', 'order']
    right_gff_df = pd.read_csv(args.rgff, sep='\t', header=None, names=right_gff_cols)
    left_gff_cols = ['chr', 'start', 'end', 'strand', 'oldid', 'id', 'order']
    left_gff_df = pd.read_csv(args.lgff, sep='\t', header=None, names=left_gff_cols)

    # 创建 GFF 类
    class GFF:
        def __init__(self, id, chr, start, end, strand, oldid, order):
            self.id = id
            self.chr = chr
            self.start = start
            self.end = end
            self.strand = strand
            self.oldid = oldid
            self.order = order

        def __repr__(self):
            return f"GFF(id={self.id}, chr={self.chr}, start={self.start}, end={self.end}, strand={self.strand}, oldid={self.oldid}, order={self.order})"

    # 将 DataFrame 转换为 GFF 对象列表
    right_gff_objects = [GFF(row['id'], row['chr'], row['start'], row['end'], row['strand'], row['oldid'], row['order'])
                         for _, row in right_gff_df.iterrows()]
    left_gff_objects = [GFF(row['id'], row['chr'], row['start'], row['end'], row['strand'], row['oldid'], row['order'])
                        for _, row in left_gff_df.iterrows()]

    # 读取 ID 文件并创建字典
    right_id_df = pd.read_csv(args.rid, sep='\t', header=None, names=['id', 'fid'])
    spRightId2fid = right_id_df.set_index('id')['fid'].to_dict()
    left_id_df = pd.read_csv(args.lid, sep='\t', header=None, names=['id', 'fid'])
    spLeftId2fid = left_id_df.set_index('id')['fid'].to_dict()

    # 验证 ID 是否存在
    assert all(id_ in right_gff_df['id'].values for id_ in spRightId2fid.keys()), "Some IDs in Right_id.txt are not present in Right.new.gff.txt"
    assert all(id_ in left_gff_df['id'].values for id_ in spLeftId2fid.keys()), "Some IDs in Left_id.txt are not present in Left.new.gff.txt"

    # 排序
    def sort_by_chr_and_end(gff_df, id_to_fid):
        data = []
        for id_ in id_to_fid.keys():
            if id_ in gff_df['id'].values:
                row = gff_df[gff_df['id'] == id_].iloc[0]
                fid = id_to_fid[id_]
                chr_ = row['chr']
                end = row['end']
                data.append((id_, fid, chr_, end))
        sorted_data = sorted(data, key=lambda x: (x[2], x[3]))
        return sorted_data

    sorted_spRight = sort_by_chr_and_end(right_gff_df, spRightId2fid)
    sorted_spLeft = sort_by_chr_and_end(left_gff_df, spLeftId2fid)

    # 聚类
    def cluster_ids(sorted_list, prefix, window_size=2e6, min_cluster_size=3):
        cluster_dict = {}
        cluster_name_dict = {}
        current_cluster = {}
        cluster_count = {}

        for i, (id_, fid, chr_, end) in enumerate(sorted_list):
            if chr_ not in cluster_count:
                cluster_count[chr_] = 1
            if chr_ not in current_cluster:
                current_cluster[chr_] = [(id_, fid, end)]
            else:
                prev_end = current_cluster[chr_][-1][2]
                if end - prev_end <= window_size:
                    current_cluster[chr_].append((id_, fid, end))
                else:
                    if len(current_cluster[chr_]) >= min_cluster_size:
                        cluster_name = f"{prefix}{chr_}_cluster_{cluster_count[chr_]}"
                        cluster_count[chr_] += 1
                    else:
                        cluster_name = "clusterUn"
                    for id__, fid__, _ in current_cluster[chr_]:
                        cluster_dict[id__] = cluster_name
                        cluster_name_dict[id__] = (fid__, cluster_name)
                    current_cluster[chr_] = [(id_, fid, end)]

        for chr_, cluster in current_cluster.items():
            if len(cluster) >= min_cluster_size:
                cluster_name = f"{prefix}{chr_}_cluster_{cluster_count[chr_]}"
            else:
                cluster_name = "clusterUn"
            for id__, fid__, _ in cluster:
                cluster_dict[id__] = cluster_name
                cluster_name_dict[id__] = (fid__, cluster_name)

        return cluster_dict, cluster_name_dict

    spRightId2cluster, spRightClusterDetails = cluster_ids(sorted_spRight, "Right_", window_size=args.window_size,
                                                           min_cluster_size=args.min_cluster_size)
    spLeftId2cluster, spLeftClusterDetails = cluster_ids(sorted_spLeft, "Left_", window_size=args.window_size,
                                                         min_cluster_size=args.min_cluster_size)

    # 输出结果
    def output_clusters(cluster_details, output_file):
        with open(output_file, 'w') as f:
            f.write("id\tfid\tcluster\n")
            for id_, (fid_, cluster) in cluster_details.items():
                f.write(f"{id_}\t{fid_}\t{cluster}\n")
        print(f"File generated: {output_file}")  # 添加提示信息

    output_clusters(spRightClusterDetails, f"{args.odir}/spRightId2cluster.txt")
    output_clusters(spLeftClusterDetails, f"{args.odir}/spLeftId2cluster.txt")


# 运行 cluster2block 功能
def run_cluster2block(args):
    if not all([args.lcluster, args.rcluster, args.bkfile]):
        print("Error: For --cluster2block, you must provide -lcluster, -rcluster, and -bkfile.")
        return

    # 加载簇文件
    spLeftId2cluster = load_cluster_file(args.lcluster)
    spRightId2cluster = load_cluster_file(args.rcluster)

    # 处理块文件
    def colinearBlock(file_path):
        with open(file_path, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
        filtered_lines = []
        block_num = None
        for line in lines:
            if line.startswith('the'):
                match = re.match(r'the (\d+)th path', line)
                if match:
                    block_num = int(match.group(1))
            elif line.startswith('overlap with block') or line.startswith('+') or line.startswith('>'):
                continue
            else:
                filtered_lines.append((line, block_num))

        data = []
        for line, block_num in filtered_lines:
            parts = line.split()
            if len(parts) >= 3:
                left_id = parts[0]
                right_id = parts[2]
                data.append([left_id, right_id, block_num])

        df = pd.DataFrame(data, columns=['leftId', 'rightId', 'blockNum'])
        return df

    block_df = colinearBlock(args.bkfile)

    # 更新数据框
    def update_dataframe_with_clusters(df, spLeftId2cluster, spRightId2cluster):
        df['leftCluster'] = df['leftId'].apply(lambda x: spLeftId2cluster.get(x, ''))
        df['rightCluster'] = df['rightId'].apply(lambda x: spRightId2cluster.get(x, ''))
        return df

    updated_block_df = update_dataframe_with_clusters(block_df, spLeftId2cluster, spRightId2cluster)

    # 保存更新后的数据框
    def save_updated_dataframe(df, output_file):
        df.to_csv(output_file, sep='\t', index=False)
        print(f"File generated: {output_file}")  # 添加提示信息

    save_updated_dataframe(updated_block_df, f"{args.odir}/left_right_block_family_cluster.txt")


# 加载簇文件
def load_cluster_file(file_path):
    cluster_dict = {}
    with open(file_path, 'r') as f:
        next(f)  # 跳过表头
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 3:
                id_, fid, cluster = parts
                cluster_dict[id_] = cluster
    return cluster_dict


if __name__ == "__main__":
    main()