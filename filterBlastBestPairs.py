import argparse

parser = argparse.ArgumentParser(description='Filter blast results and write to an output file.')
parser.add_argument('-bl', '--blast_file', type=str, help='Path to the blast file', required=True)
parser.add_argument('-o', '--output_file', type=str, help='Path to the output file', required=True)
args = parser.parse_args()

best_hit_num = 1

try:
    with open(args.blast_file, 'r') as bl:
        with open(args.output_file, 'w') as out:
            subject_num = 0
            query = ""
            last_query = ""
            last_subject = ""
            hit_num = best_hit_num

            for line in bl:
                array = line.split()
                if len(array) < 2:
                    continue
                if array[0] == last_query and array[1] == last_subject:
                    continue
                if array[0] == last_query:
                    if hit_num >= best_hit_num:
                        continue
                    hit_num += 1

                out.write(array[0] + "\t" + array[1] + "\n")
                last_query = array[0]
                last_subject = array[1]

except IOError as e:
    print("I/O error:", e)