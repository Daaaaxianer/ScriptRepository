# -*- coding: utf-8 -*-
# FileName  : extend_wgd_scope.py
# Time      : 2025/4/26 20:26
# Author    : xian

import sys
import re
import os

def read_wgd_ids(wgd_file):
    """
    读取WGD文件中的所有基因ID（所有列）
    参数:
        wgd_file: WGD文件路径
    返回:
        包含所有WGD基因ID的集合（所有列合并后去重）
    """
    wgd_ids = set()
    with open(wgd_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            for gene_id in parts:
                if gene_id:
                    wgd_ids.add(gene_id)
    return wgd_ids

def parse_gff_file(gff_file):
    """
    解析GFF文件并建立索引结构
    参数:
        gff_file: GFF文件路径, 默认格式：chr, start, end, strand, oldid, id, order
    返回:
        tuple: (id_to_info, chr_order_map)
        id_to_info: 字典，基因ID到(chr, order)的映射
        chr_order_map: 字典，染色体到order到基因ID列表的映射
    """
    id_to_info = {}
    chr_order_map = {}
    
    with open(gff_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            
            if len(parts) < 7:
                continue
                
            chr_num = parts[0]
            gene_id = parts[5]
            
            try:
                order = int(parts[6])
            except ValueError:
                continue
            
            id_to_info[gene_id] = (chr_num, order)
            
            if chr_num not in chr_order_map:
                chr_order_map[chr_num] = {}
            if order not in chr_order_map[chr_num]:
                chr_order_map[chr_num][order] = []
            chr_order_map[chr_num][order].append(gene_id)
    
    return id_to_info, chr_order_map

def find_adjacent_genes(wgd_ids, id_to_info, chr_order_map):
    """
    查找每个WGD基因前后50个order范围内的所有基因
    参数:
        wgd_ids: 需要查找的WGD基因ID集合
        id_to_info: 基因ID到(chr, order)的映射
        chr_order_map: 染色体到order到基因ID列表的映射
    返回:
        集合: 包含所有符合条件的基因ID（已去重）
        字典: 基因ID到(chr, order)的映射（用于排序）
    """
    result_ids = set()
    found_ids = set()
    not_found_ids = set()
    id_to_chr_order = {}  # 用于存储结果基因的chr和order信息
    
    for wgd_id in wgd_ids:
        if wgd_id not in id_to_info:
            not_found_ids.add(wgd_id)
            continue
            
        chr_num, order = id_to_info[wgd_id]
        found_ids.add(wgd_id)
        
        if chr_num not in chr_order_map:
            continue
        
        chr_orders = sorted(chr_order_map[chr_num].keys())
        
        try:
            current_idx = chr_orders.index(order)
        except ValueError:
            continue
        
        start_idx = max(0, current_idx - 50)
        end_idx = min(len(chr_orders) - 1, current_idx + 50)
        
        for o in chr_orders[start_idx:end_idx + 1]:
            if o in chr_order_map[chr_num]:
                for gene_id in chr_order_map[chr_num][o]:
                    if gene_id not in result_ids:  # 去重
                        result_ids.add(gene_id)
                        id_to_chr_order[gene_id] = (chr_num, o)
    
    print(f"成功匹配 {len(found_ids)} 个WGD ID")
    if not_found_ids:
        print(f"未找到 {len(not_found_ids)} 个WGD ID，示例: {', '.join(list(not_found_ids)[:5])}{'...' if len(not_found_ids)>5 else ''}")
    
    return result_ids, id_to_chr_order

def natural_sort_key(s):
    """
    自然排序键函数，用于chr的排序（如chr1, chr2,...chr10）
    """
    return [int(text) if text.isdigit() else text.lower() 
            for text in re.split('([0-9]+)', s)]

def write_sorted_output(output_file, gene_ids, id_to_chr_order):
    """
    将结果按chr和order排序后写入输出文件
    参数:
        output_file: 输出文件路径
        gene_ids: 基因ID集合
        id_to_chr_order: 基因ID到(chr, order)的映射
    """
    # 转换为列表并按chr和order排序
    sorted_genes = sorted(
        gene_ids,
        key=lambda x: (natural_sort_key(id_to_chr_order[x][0]), id_to_chr_order[x][1])
    )
    
    with open(output_file, 'w') as f:
        for gene_id in sorted_genes:
            f.write(f"{gene_id}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        script_name = os.path.basename(sys.argv[0])
        print(f"\n用法: python {script_name} <gff文件> <wgd文件> <输出文件>\n")
        print(f"示例: python {script_name} xxx.gff xxx.WGD.txt xxx.WGD.extend.txt\n")
        sys.exit(1)
    
    gff_file = sys.argv[1]
    wgd_file = sys.argv[2]
    output_file = sys.argv[3]
    
    print("正在读取WGD基因ID...")
    wgd_ids = read_wgd_ids(wgd_file)
    print(f"从WGD文件中读取到 {len(wgd_ids)} 个唯一基因ID")
    
    print("正在解析GFF文件并建立索引...")
    id_to_info, chr_order_map = parse_gff_file(gff_file)
    print(f"GFF文件解析完成，共找到 {len(id_to_info)} 个基因")
    
    print("正在查找相邻基因并去重...")
    result_ids, id_to_chr_order = find_adjacent_genes(wgd_ids, id_to_info, chr_order_map)
    print(f"找到 {len(result_ids)} 个唯一符合条件的基因ID")
    
    print("正在按chr和order排序并写入输出文件...")
    write_sorted_output(output_file, result_ids, id_to_chr_order)
    print(f"处理完成，结果已保存到 {output_file}")

# if __name__ == "__main__":
#     main()