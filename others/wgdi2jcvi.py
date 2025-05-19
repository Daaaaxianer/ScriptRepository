# -*- coding:UTF-8 -*-
# FileName  :wgdi2jcvi.py
# Time      :2025/1/09
# Author    :xian


import argparse
import pandas as pd


def load_filter_ids(filter_filepath):
    with open(filter_filepath, 'r') as file:
        filter_ids = set(int(line.strip()) for line in file if line.strip())
    return filter_ids


class Alignment:
    def __init__(self, head_info, alignment_id, score, pvalue, length, chrpair, direction):
        self.head = head_info
        self.id = alignment_id
        self.score = score
        self.pvalue = pvalue
        self.length = length
        self.chrpair = chrpair
        self.dirct = direction
        self.bkpair = []

    def add_bkpair(self, bkpair):
        self.bkpair.append(bkpair)

    def __repr__(self):
        return (f"Alignment(id={self.id}, head='{self.head}', score={self.score}, "
                f"pvalue={self.pvalue}, length={self.length}, chrpair={self.chrpair}, dirct={self.dirct})")


def parse_file(filepath, filter_filepath=None):
    alignments = []
    filter_ids = load_filter_ids(filter_filepath) if filter_filepath else None

    with open(filepath, 'r') as file:
        lines = file.readlines()
        alignment = None
        for line in lines:
            if line.startswith("# Alignment"):
                parts = line.split()[2:]
                id = int(parts[0].replace(':', ''))

                if not filter_ids or id in filter_ids:
                    if alignment:
                        alignments.append(alignment)
                    score = int(parts[1].split('=')[1])
                    pvalue = float(parts[2].split('=')[1])
                    length = int(parts[3].split('=')[1])
                    chrpair = parts[4]
                    direction = '+' if parts[5].lower() == 'plus' else '-'  # 在这里进行转换
                    alignment = Alignment(line.strip(), id, score, pvalue, length, chrpair, direction)
                else:
                    alignment = None
            elif line.strip() and not line.startswith('#'):
                if alignment:
                    bkpair_parts = line.split()
                    bkpair_info = {
                        "query_id": bkpair_parts[0],
                        "query_num": int(bkpair_parts[1]),
                        "subject_id": bkpair_parts[2],
                        "subject_num": int(bkpair_parts[3]),
                        "dir": int(bkpair_parts[4])
                    }
                    alignment.add_bkpair(bkpair_info)
        if alignment:
            alignments.append(alignment)
    return alignments


def save_alignment_summary_with_pandas(alignments, output_filepath):
    data = []
    for alignment in alignments:
        first_query_id = alignment.bkpair[0]['query_id'] if alignment.bkpair else 'N/A'
        last_query_id = alignment.bkpair[-1]['query_id'] if alignment.bkpair else 'N/A'
        first_subject_id = alignment.bkpair[0]['subject_id'] if alignment.bkpair else 'N/A'
        last_subject_id = alignment.bkpair[-1]['subject_id'] if alignment.bkpair else 'N/A'

        data.append([first_query_id, last_query_id, first_subject_id, last_subject_id,
                     alignment.length, alignment.dirct])

    df = pd.DataFrame(data)

    df.to_csv(output_filepath, sep='\t', header=False, index=False)

def save_alignment_pairs(alignments, output_filepath):
    data = []
    for alignment in alignments:
        for bkpair in alignment.bkpair:
            data.append([bkpair['query_id'], bkpair['subject_id']])

    df = pd.DataFrame(data)

    df.to_csv(output_filepath, sep='\t', header=False, index=False)

def gff2bed(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t', header=None)

    df_transformed = df.loc[:, [0, 2, 3, 1, 4, 6]]
    df_transformed.insert(4, 'fixed_value', '0')

    cols = [0, 2, 3, 1, 'fixed_value', 4]
    df_transformed = df_transformed[cols]

    df_transformed.to_csv(output_file, sep='\t', header=False, index=False)


def load_blast_results(blast_file):
    blast_dict = {}
    with open(blast_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            query = parts[0]
            subject = parts[1]
            score = int(float(parts[-1]))  # 将得分转换为整数

            # 确保每个query-subject对只有一个得分
            if query not in blast_dict:
                blast_dict[query] = {}
            blast_dict[query][subject] = score

    return blast_dict


def generate_anchor_file(alignments, blast_dict, output_file):
    anchor_data = []
    for alignment in alignments:
        for bkpair in alignment.bkpair:
            query_id = bkpair['query_id']
            subject_id = bkpair['subject_id']

            # 获取query和subject的得分并转换为整数
            query_score = blast_dict.get(query_id, {}).get(subject_id, 'N/A')
            if query_score != 'N/A':
                query_score = int(query_score)

            anchor_data.append([query_id, subject_id, query_score])

    df = pd.DataFrame(anchor_data)

    df.to_csv(output_file, sep='\t', header=False, index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Process collinearity file, convert GFF to BED, or generate anchors from BLAST results.")
    subparsers = parser.add_subparsers(dest='command')

    # --blast2anchor
    parser_blast2anchor = subparsers.add_parser('blast2anchor',
                                                help='Generate anchor file based on collinearity and BLAST results.')
    parser_blast2anchor.add_argument('-c', required=True, metavar='collinearity_file',
                                     help='Path to the input collinearity file.')
    parser_blast2anchor.add_argument('-b', required=True, metavar='blast_file', help='Path to the BLAST result file.')
    parser_blast2anchor.add_argument('-o', required=True, metavar='output_anchor_file',
                                     help='Path to the output anchor file.')

    # --gff2bed
    parser_gff2bed = subparsers.add_parser('gff2bed', help='Convert GFF file to BED format.')
    parser_gff2bed.add_argument('-gff', required=True, metavar='input_gff_file', help='Path to the input GFF file.')
    parser_gff2bed.add_argument('-bed', required=True, metavar='output_bed_file', help='Path to the output BED file.')

    # --block2simple
    parser_block2simple = subparsers.add_parser('block2simple',
                                                help='Process collinearity file and output a simple summary.')
    parser_block2simple.add_argument('-c', required=True, metavar='input_collinearity_file',
                                     help='Path to the input collinearity file.')
    parser_block2simple.add_argument('-o', required=True, metavar='output_simple_file',
                                     help='Path to the output simple file.')
    parser_block2simple.add_argument('--filter', metavar='filter_blockid_file',
                                     help='Path to the filter block ID file.')

    # --gff2pair
    parser_block2pair = subparsers.add_parser('block2pair', help='Extract gene pairs of certain blocks from a collinearity file.')
    parser_block2pair.add_argument('-c', required=True, metavar='input_collinearity_file',
                                     help='Path to the input collinearity file.')
    parser_block2pair.add_argument('--filter', metavar='filter_blockid_file',
                                     help='Path to the filter block ID file.')
    parser_block2pair.add_argument('-o', required=True, metavar='output_pair_file',
                                     help='Path to the output pair file.')

    args = parser.parse_args()

    if args.command == 'blast2anchor':
        alignments = parse_file(args.c)
        blast_dict = load_blast_results(args.b)
        generate_anchor_file(alignments, blast_dict, args.o)
    elif args.command == 'gff2bed':
        gff2bed(args.gff, args.bed)
    elif args.command == 'block2simple':
        filter_filepath = args.filter if args.filter else None
        alignments = parse_file(args.c, filter_filepath)
        save_alignment_summary_with_pandas(alignments, args.o)
    elif args.command == 'block2pair':
        filter_filepath = args.filter if args.filter else None
        alignments = parse_file(args.c, filter_filepath)
        save_alignment_pairs(alignments, args.o)
    else:
        print("Please specify one of the commands: blast2anchor, gff2bed, block2simple.")


if __name__ == "__main__":
    main()