# -*- coding: utf-8 -*-
# FileName  : wgdi_bkinfo_to_circos_link.py
# Time      : 2025/1/14 13:09
# Author    : xian

import pandas as pd


def create_index(gff_file):
    """
    从gff文件中读取数据并创建一个多维字典。
    :param gff_file: 包含染色体、ID、起始位置、终止位置、链方向、顺序和旧ID的文件路径。
    :return: 多维字典，可以通过chr索引order，并进一步通过order索引start和end。
    """
    # 读取gff文件
    df = pd.read_csv(gff_file, sep='\t', header=None, names=['chr', 'id', 'start', 'end', 'strand', 'order', 'oldid'])
    # df = pd.read_csv(gff_file, sep='\t', header=None, names=['chr', 'start', 'end', 'strand', 'oldid', 'id','order'])

    index = {}
    for _, row in df.iterrows():
        if row['chr'] not in index:
            index[row['chr']] = {}
        index[row['chr']][row['order']] = {'start': row['start'], 'end': row['end']}
    return index


def replace_values(chr_key, pos_key, order, index):
    """
    根据给定的染色体编号（chr_key）、位置类型（pos_key）和顺序编号（order），从index中查找对应的起始或终止位置。
    :param chr_key: 染色体编号。
    :param pos_key: 位置类型，可以是'start'或'end'。
    :param order: 顺序编号。
    :param index: 多维字典。
    :return: 对应的起始或终止位置，如果没有找到匹配项，则返回None。
    """
    if chr_key in index and order in index[chr_key]:
        return index[chr_key][order][pos_key]
    else:
        print(f"Warning: No match found for chr: {chr_key}, order: {order}")
        return None


def process_block_info(block_info_file, index, output_file):
    """
    处理wgdi生成的blockinfo文件，利用index字典替换start1, end1, start2, end2字段，并将结果写入新的文件。
    :param block_info_file: 需要处理的输入文件路径。
    :param index: 多维字典，用于查找对应的起始和终止位置。
    :param output_file: 输出文件路径。
    """
    # 读取blockinfo文件
    df = pd.read_csv(block_info_file, sep=',', header=0)
    print("All columns:", df.columns.tolist())
    
    # 明确地创建df_block作为一个独立的副本
    df_block = df[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']].copy()
    print(df_block)

    # 使用.loc进行赋值操作
    df_block.loc[:, 'start1'] = df_block.apply(lambda x: replace_values(x['chr1'], 'start', x['start1'], index), axis=1)
    df_block.loc[:, 'end1'] = df_block.apply(lambda x: replace_values(x['chr1'], 'end', x['end1'], index), axis=1)
    df_block.loc[:, 'start2'] = df_block.apply(lambda x: replace_values(x['chr2'], 'start', x['start2'], index), axis=1)
    df_block.loc[:, 'end2'] = df_block.apply(lambda x: replace_values(x['chr2'], 'end', x['end2'], index), axis=1)

    # 将修改后的数据框写入新的文件，并使用制表符作为分隔符
    df_block.to_csv(output_file, sep='\t', header=False, index=False)


if __name__ == "__main__":
    import sys

    # 检查命令行参数数量
    if len(sys.argv) != 4:
        print("Usage: python script.py coca.new.gff.txt coca_coca.blockinfo.notandem.txt output.txt")
        sys.exit(1)

    # 获取命令行参数
    gff_file = sys.argv[1]  # gff文件路径
    block_info_file = sys.argv[2]  # blockinfo文件路径
    output_file = sys.argv[3]  # 输出文件路径

    # 创建索引字典
    index = create_index(gff_file)

    # 处理blockinfo文件并生成输出文件
    process_block_info(block_info_file, index, output_file)