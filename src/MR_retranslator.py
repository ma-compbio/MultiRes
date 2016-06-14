import sys


def translator(car_content,gene_dict):
    return [gene_dict[x] for x in car_content]


gene_id_dict = {}

with open(sys.argv[1]) as gene_id_file:
    for line in gene_id_file:
        if not(line.isspace()):
            g, id = line[1:].strip().split(' ')
            gene_id_dict[id] = g


with open(sys.argv[2]) as pq_tree_file:
    for line in pq_tree_file:
        if 'Q' not in line:
            print line.strip()
        else:
            content = line.strip().split(' ')[1:-1]
            translated = translator(content,gene_id_dict)
            print '_Q '+' '.join(translated)+' Q_'
