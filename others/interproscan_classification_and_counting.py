#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : yuzijian
# @Email         : yuzijian1010@163.com
# @FileName      : interproscan_gene_function_count.py
# @Time          : 2024-09-04 22:21:58
# @description   : interproscan 注释出来的功能基因不同库的数目统计
"""

import argparse


def process_tsv_file(input_file, output_file):
    column_4_counts = {}
    gene_annotations = set()
    processed_pairs = set()

    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) > 3:
                gene_id = parts[0]
                column_4_value = parts[3]
                pair = (gene_id, column_4_value)
                if pair in processed_pairs:
                    continue
                if column_4_value in column_4_counts:
                    column_4_counts[column_4_value] += 1
                else:
                    column_4_counts[column_4_value] = 1
                processed_pairs.add(pair)
                gene_annotations.add(gene_id)

    with open(output_file, 'w') as output:
        for key, value in column_4_counts.items():
            output.write(f"{key}\t{value}\n")

    print(f"第一列去重后的基因注释个数为: {len(gene_annotations)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process interproscan result file and count gene function numbers.')
    parser.add_argument('input_file', type=str, help='Path to the input interproscan result file (in tsv format)')
    parser.add_argument('output_file', type=str, help='Path to the output file for storing the count results')
    args = parser.parse_args()
    process_tsv_file(args.input_file, args.output_file)
    