#!/usr/bin/env

import sys
import DataStructures.markers as markers
import DataStructures.genomes as genomes

def adjacentPos(gen,chrom,st_pos):
    i = gen.findLocus(chrom,st_pos)
    if i%2 == 0:
        if i == 0:
            return ''
        return gen.chromosomes[chrom][i-1].name()
    else:
        if i == len(gen.chromosomes[chrom]) or i == len(gen.chromosomes[chrom])-1:
            return ''
        return gen.chromosomes[chrom][i+1].name()

def inSpecies(families,sp_name,gen,m0,m1):
    for loc in families[m0].loci:
        if loc.species == sp_name:
            chrom = loc.chromosome
            st_pos = loc.start
            end_pos = loc.end
            if adjacentPos(gen,chrom,st_pos) == m1:
                return True
    return False

def conservedAdjs(gens,families,pair,adj_dict):
    ref_sp = gens[pair[0]].species
    compare_sp = gens[pair[1]].species
    ref_gen = gens[pair[0]]
    other_gen = gens[pair[1]]
    for chrom in ref_gen.chromosomes:
        current_chrom = ref_gen.chromosomes[chrom]
        for i in range(1,len(current_chrom)-1,2):
            m0 = current_chrom[i].name()
            m1 = current_chrom[i+1].name()
            adj_pair = tuple(sorted([m0,m1]))
            if inSpecies(families,compare_sp,other_gen,m0,m1):
                if adj_pair not in adj_dict:
                    adj_dict[adj_pair] = []
                if ref_sp not in adj_dict[adj_pair]:
                    adj_dict[adj_pair].append(ref_sp)
                if compare_sp not in adj_dict[adj_pair]:
                    adj_dict[adj_pair].append(compare_sp)
    return adj_dict

def getAdjacencies(all_genomes,families,sp_pairs):
    adj_dictionary = {} # Dictionary of adjacency tuples, mapping to species list containing them
    for pair in sp_pairs:
        adj_dictionary.update(conservedAdjs(all_genomes,families,pair,adj_dictionary))
    return adj_dictionary


def read_pairs_file(filename):
    pairs = []
    with open(filename) as pairs_file:
        for line in pairs_file:
            if not(line.isspace()) and '#' not in line:
                pair = sorted(line.strip().split(' '))
                pairs.append(pair)
    return sorted(pairs)

if __name__ == "__main__" :
    families = markers.read_hom_families_file(sys.argv[1])
    sp_pairs = read_pairs_file(sys.argv[2])
    all_sp = []
    for x in sp_pairs:
        all_sp += x
    unique_sp = list(set(all_sp))
    extant_genomes = genomes.getGenomes(families,unique_sp)
    adjacencies = getAdjacencies(extant_genomes,families,sp_pairs)
    i = 0
    for a in adjacencies:
        print str(i)+"|1;"+','.join(sorted(adjacencies[a]))+':'+' '.join(a)
        i += 1
