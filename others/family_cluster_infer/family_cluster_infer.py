# -*- coding:UTF-8 -*-
# FileName  :family_cluster_infer.py
# Time      :2025/2/10 
# Author    :xian


import pandas as pd
import re
import argparse

import pandas as pd
import re
import argparse


# 添加命令行参数解析
def parse_args():
    parser = argparse.ArgumentParser(description="Process GFF and block files for clustering.")
    parser.add_argument("--right_gff", type=str, required=True, help="Path to Right.new.gff.txt file")
    parser.add_argument("--left_gff", type=str, required=True, help="Path to Left.new.gff.txt file")
    parser.add_argument("--right_id", type=str, required=True, help="Path to Right_id.txt file")
    parser.add_argument("--left_id", type=str, required=True, help="Path to Left_id.txt file")
    parser.add_argument("--block_file", type=str, required=True, help="Path to Left_Right.block.rr.txt file")
    parser.add_argument("--window_size", type=float, default=2e6,
                        help="Sliding window size in base pairs (default: 2e6)")
    parser.add_argument("--min_cluster_size", type=int, default=3, help="Minimum cluster size (default: 3)")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory path for output files")
    return parser.parse_args()


# 主函数
def main():
    # 解析命令行参数
    args = parse_args()

    print("Starting script with the following parameters:")
    print(f"Right GFF File: {args.right_gff}")
    print(f"Left GFF File: {args.left_gff}")
    print(f"Right ID File: {args.right_id}")
    print(f"Left ID File: {args.left_id}")
    print(f"Block File: {args.block_file}")
    print(f"Output Directory: {args.output_dir}")
    print(f"Sliding Window Size: {args.window_size} bp")
    print(f"Minimum Cluster Size: {args.min_cluster_size}")

    # 1. 读取 Right.new.gff.txt 文件
    right_gff_cols = ['chr', 'start', 'end', 'strand', 'oldid', 'id', 'order']
    right_gff_df = pd.read_csv(args.right_gff, sep='\t', header=None, names=right_gff_cols)
    print("Right GFF file loaded successfully.")

    # 2. 读取 Left.new.gff.txt 文件
    left_gff_cols = ['chr', 'start', 'end', 'strand', 'oldid', 'id', 'order']
    left_gff_df = pd.read_csv(args.left_gff, sep='\t', header=None, names=left_gff_cols)
    print("Left GFF file loaded successfully.")

    # 3. 创建 GFF 类
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
    print("GFF objects created successfully.")

    # 4. 读取 Right_id.txt 并创建字典 spRightId2fid
    right_id_df = pd.read_csv(args.right_id, sep='\t', header=None, names=['id', 'fid'])
    spRightId2fid = right_id_df.set_index('id')['fid'].to_dict()
    print("Right ID file processed successfully.")

    # 5. 读取 Left_id.txt 并创建字典 spLeftId2fid
    left_id_df = pd.read_csv(args.left_id, sep='\t', header=None, names=['id', 'fid'])
    spLeftId2fid = left_id_df.set_index('id')['fid'].to_dict()
    print("Left ID file processed successfully.")

    # 6. 确保文件中的 id 存在于对应的 .gff 文件中
    assert all(id_ in right_gff_df['id'].values for id_ in
               spRightId2fid.keys()), "Some IDs in Right_id.txt are not present in Right.new.gff.txt"
    assert all(id_ in left_gff_df['id'].values for id_ in
               spLeftId2fid.keys()), "Some IDs in Left_id.txt are not present in Left.new.gff.txt"
    print("ID validation completed successfully.")

    # 7. 按照 id.chr 和 id.end 排序
    def sort_by_chr_and_end(gff_objects, id_to_fid):
        sorted_list = []
        chr_groups = {}

        # 将 id 对应的 fid 添加到 GFF 对象中
        for obj in gff_objects:
            if obj.id in id_to_fid:
                fid = id_to_fid[obj.id]
                if obj.chr not in chr_groups:
                    chr_groups[obj.chr] = []
                chr_groups[obj.chr].append((obj, fid))

        # 对每个染色体组按 end 排序
        for chr_, group in chr_groups.items():
            sorted_group = sorted(group, key=lambda x: x[0].end)
            sorted_list.extend(sorted_group)

        return sorted_list

    sorted_spRight = sort_by_chr_and_end(right_gff_objects, spRightId2fid)
    sorted_spLeft = sort_by_chr_and_end(left_gff_objects, spLeftId2fid)
    print("Sorting by chromosome and end position completed successfully.")

    # 8. 使用滑动窗口进行 cluster 划分
    def cluster_ids(sorted_list, prefix, window_size=2e6, min_cluster_size=3):
        cluster_dict = {}
        cluster_name_dict = {}
        current_cluster = {}
        cluster_count = {}

        for i, (gff_obj, fid) in enumerate(sorted_list):
            chr_ = gff_obj.chr
            end = gff_obj.end

            if chr_ not in cluster_count:
                cluster_count[chr_] = 1

            if chr_ not in current_cluster:
                current_cluster[chr_] = [(gff_obj.id, fid, end)]
            else:
                prev_end = current_cluster[chr_][-1][2]
                if end - prev_end <= window_size:  # 使用滑动窗口大小作为参数
                    current_cluster[chr_].append((gff_obj.id, fid, end))
                else:
                    # 检查当前 cluster 是否满足最小 size
                    if len(current_cluster[chr_]) >= min_cluster_size:
                        cluster_name = f"{prefix}{chr_}_cluster_{cluster_count[chr_]}"
                        cluster_count[chr_] += 1
                    else:
                        cluster_name = "clusterUn"

                    # 更新字典
                    for id_, fid_, _ in current_cluster[chr_]:
                        cluster_dict[id_] = cluster_name
                        cluster_name_dict[id_] = (fid_, cluster_name)

                    # 初始化新的 cluster
                    current_cluster[chr_] = [(gff_obj.id, fid, end)]

        # 处理最后一个 cluster
        if current_cluster:
            chr_ = list(current_cluster.keys())[0]
            if len(current_cluster[chr_]) >= min_cluster_size:
                cluster_name = f"{prefix}{chr_}_cluster_{cluster_count[chr_]}"
            else:
                cluster_name = "clusterUn"

            for id_, fid_, _ in current_cluster[chr_]:
                cluster_dict[id_] = cluster_name
                cluster_name_dict[id_] = (fid_, cluster_name)

        return cluster_dict, cluster_name_dict

    # 调用 cluster_ids 函数
    spRightId2cluster, spRightClusterDetails = cluster_ids(sorted_spRight, "Right_", window_size=args.window_size,
                                                           min_cluster_size=args.min_cluster_size)
    spLeftId2cluster, spLeftClusterDetails = cluster_ids(sorted_spLeft, "Left_", window_size=args.window_size,
                                                         min_cluster_size=args.min_cluster_size)
    print("Clustering completed successfully.")

    # 输出结果到文件
    def output_clusters(cluster_details, output_file):
        with open(output_file, 'w') as f:
            f.write("id\tfid\tcluster\n")
            for id_, (fid_, cluster) in cluster_details.items():
                f.write(f"{id_}\t{fid_}\t{cluster}\n")

    output_clusters(spRightClusterDetails, f"{args.output_dir}/spRightId2cluster.txt")
    output_clusters(spLeftClusterDetails, f"{args.output_dir}/spLeftId2cluster.txt")
    print("Cluster results saved to spRightId2cluster.txt and spLeftId2cluster.txt.")

    # 9. 定义 colinearBlock 函数
    def colinearBlock(file_path):
        # 读取文件并过滤行
        with open(file_path, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]

        filtered_lines = []
        block_num = None

        for line in lines:
            if line.startswith('the'):
                # 使用正则表达式提取 the 前的数字
                match = re.match(r'the (\d+)th path', line)
                if match:
                    block_num = int(match.group(1))  # 提取数字
            elif line.startswith('overlap with block'):
                # 跳过以 "overlap with block" 开头的行
                continue
            elif line.startswith('+') or line.startswith('>'):
                # 跳过以 "+" 或 ">" 开头的行
                continue
            else:
                # 保留其他行，并将当前 block_num 与该行关联
                filtered_lines.append((line, block_num))

        # 创建 DataFrame
        data = []
        for line, block_num in filtered_lines:
            parts = line.split()
            if len(parts) >= 3:
                left_id = parts[0]
                right_id = parts[2]
                data.append([left_id, right_id, block_num])

        df = pd.DataFrame(data, columns=['leftId', 'rightId', 'blockNum'])
        return df

    # 10. 更新新数据框
    def update_dataframe_with_clusters(df, spLeftId2cluster, spRightId2cluster):
        # 添加 leftCluster 和 rightCluster 列
        df['leftCluster'] = df['leftId'].apply(lambda x: spLeftId2cluster.get(x, ''))
        df['rightCluster'] = df['rightId'].apply(lambda x: spRightId2cluster.get(x, ''))
        return df

    # 11. 输出更新后的数据框到新文件
    def save_updated_dataframe(df, output_file):
        df.to_csv(output_file, sep='\t', index=False)

    # 调用 colinearBlock 函数处理 Left_Right.block.rr.txt 文件
    block_df = colinearBlock(args.block_file)
    print("Block file processed successfully.")

    # 更新数据框
    updated_block_df = update_dataframe_with_clusters(block_df, spLeftId2cluster, spRightId2cluster)
    print("Dataframe updated with clusters successfully.")

    # 保存更新后的数据框到文件
    save_updated_dataframe(updated_block_df, f"{args.output_dir}/left_right_block_family_cluster.txt")
    print("Updated block data saved to left_right_block_family_cluster.txt.")
    print("Script execution completed successfully.")


if __name__ == "__main__":
    main()
