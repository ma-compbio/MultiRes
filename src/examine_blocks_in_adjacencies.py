#!/usr/bin/env
import sys
from marker_string import *
from DataStructures import genomes
from DataStructures import markers

species_list_in = []
species_list_out = []
with open(sys.argv[3]) as sp_list_file:
    toggle_type = False
    for line in sp_list_file:
        if "INGROUPS" in line.upper():
            toggle_type = True
        elif "OUTGROUPS" in line.lower():
            toggle_type = False
        elif not(line.isspace()):
            if toggle_type:
                if not ('#' in line):
                    species_list_in.append(line.strip())
            else:
                if not ('#' in line):
                    species_list_out.append(line.strip())

synteny_blocks_big = markers.read_hom_families_file(sys.argv[1])
synteny_blocks_small = markers.read_hom_families_file(sys.argv[2])
#uu_synteny_blocks_100k = markers.read_hom_families_file("../Mammalian_data/euth.100K.Conserved.Segments.Species_filtered.uu")
#uu_synteny_blocks_300k = markers.read_hom_families_file("../Mammalian_data/euth.300K.Conserved.Segments.Species_filtered.uu")
#uu_synteny_blocks_500k = markers.read_hom_families_file("../Mammalian_data/euth.500K.Conserved.Segments.Species_filtered.uu")
#uu_synteny_blocks_20k = markers.read_hom_families_file("../Mammalian_data/Conserved.Segments.20K.Species_filtered.uu")
#uu_synteny_blocks_50k = markers.read_hom_families_file("../Mammalian_data/Conserved.Segments.50K.Species_filtered.uu")

#synteny_blocks_100k = markers.read_hom_families_file("../Mammalian_data/euth.100K.Conserved.Segments.Species_filtered")
#synteny_blocks_300k = markers.read_hom_families_file("../Mammalian_data/euth.300K.Conserved.Segments.Species_filtered")
#synteny_blocks_500k = markers.read_hom_families_file("../Mammalian_data/euth.500K.Conserved.Segments.Species_filtered")
#synteny_blocks_20k = markers.read_hom_families_file("../Mammalian_data/Conserved.Segments.20K.Species_filtered")
#synteny_blocks_50k = markers.read_hom_families_file("../Mammalian_data/Conserved.Segments.50K.Species_filtered")

#synteny_blocks_big = {}
#
#for m in synteny_blocks_list_big:
#    print m
#    synteny_blocks_big[m.id] = m
#
#synteny_blocks_small = {}
#
#for m in synteny_blocks_list_small:
#    synteny_blocks_small[m.id] = m
#

all_species = species_list_in+species_list_out

#uu_genomes_100k = genomes.getGenomes(uu_synteny_blocks_100k,all_species)
#uu_genomes_300k = genomes.getGenomes(uu_synteny_blocks_300k,all_species)
#uu_genomes_500k = genomes.getGenomes(uu_synteny_blocks_500k,all_species)
#uu_genomes_50k = genomes.getGenomes(uu_synteny_blocks_50k,all_species)
#uu_genomes_20k = genomes.getGenomes(uu_synteny_blocks_20k,all_species)

#genomes_100k = genomes.getGenomes(synteny_blocks_100k,all_species)
#genomes_300k = genomes.getGenomes(synteny_blocks_300k,all_species)
#genomes_500k = genomes.getGenomes(synteny_blocks_500k,all_species)
#genomes_20k = genomes.getGenomes(synteny_blocks_20k,all_species)
#genomes_50k = genomes.getGenomes(synteny_blocks_50k,all_species)

genomes_big = genomes.getGenomes(synteny_blocks_big,all_species)
genomes_small  = genomes.getGenomes(synteny_blocks_small,all_species)

def getContainment(genome,chromosome,start,end):
    if start > end:
        start, end = end,start
    index = genome.findLocusFuzzy(chromosome,start)
    temp_list = []
    first_marker = True
    while index != -sys.maxint+1 and index < len(genome.chromosomes[chromosome]) and genome.chromosomes[chromosome][index].start() < end:
        orient = ''
        if end-start+1 < genome.chromosomes[chromosome][index].end()-genome.chromosomes[chromosome][index].start()+1:
            return [] 
        if genome.chromosomes[chromosome][index].orientation() < 0:
            orient = '-'
        if first_marker and abs(genome.chromosomes[chromosome][index].end()-start) >= 0.5*(genome.chromosomes[chromosome][index].end()-genome.chromosomes[chromosome][index].start()):
            temp_list.append(orient+genome.chromosomes[chromosome][index].name())
            first_marker = False
        elif first_marker and not(abs(genome.chromosomes[chromosome][index].end()-start) >= 0.5*(genome.chromosomes[chromosome][index].end()-genome.chromosomes[chromosome][index].start())):
            first_marker = False
        elif not(first_marker) and genome.chromosomes[chromosome][index].end() <= end:
            temp_list.append(orient+genome.chromosomes[chromosome][index].name())
        elif not(first_marker) and genome.chromosomes[chromosome][index].end() > end:
            if abs(end-genome.chromosomes[chromosome][index].start()) > 0.5*(genome.chromosomes[chromosome][index].end()-genome.chromosomes[chromosome][index].start()):
                temp_list.append(orient+genome.chromosomes[chromosome][index].name())
        index += 1
    return temp_list

c1p_adjacencies = []
intercar_pos = -500000
sp_list = []
with open(sys.argv[4]) as adj_file:
    index = 0#int(line[:line.find('|')])
    for line in adj_file:
        if not(line.isspace() or line[0] == "#"):
            l = line[line.find(":")+1:line.find('#')].strip().split(' ')
            sp = line[line.find(";")+1:line.find(':')].split(',')
            c1p_adjacencies.append(tuple(l))
            to_remove = []
            for s in sp:
                if s not in all_species:
                    to_remove.append(s)
            for s in to_remove:
                sp.remove(s)
            sp_list.append(sp)
#            print l
#            sys.exit(0)
#            a1 = l[0].split('_')
#            if a1[0][0] == '-' or a1[0][0] == '+':
#                a1[0] = a1[0][1:]
#            a2 = l[1].split('_')
#            if a2[0][0] == '-' or a2[0][0] == '+':
#                a2[0] = a2[0][1:]
#            a1 = tuple(a1)
#            a2 = tuple(a2)
#            if (a2,a1) not in c1p_adjacencies:
#               c1p_adjacencies.append((a1,a2))
#               index += 1
#        if "INTERCAR" in line:
#            intercar_pos = index

def findGapPositions(genome,m1,orient1,m2,orient2,first_locus):
    chr = first_locus.chromosome
    orientation = first_locus.orientation
    index = genome.findLocus(chr,first_locus.start)
    if index != -sys.maxint+1:
        if orientation > 0:
            if orient1 == 'h' and index < len(genome.chromosomes[chr])-1:
                if genome.chromosomes[chr][index+1].name() == m2:
                    if orient2 == 't' and genome.chromosomes[chr][index+1].orientation() > 0:
                        return genome.chromosomes[chr][index].end(), genome.chromosomes[chr][index+1].start(),'+'
                    elif orient2 == 'h' and genome.chromosomes[chr][index+1].orientation() < 0 :
                        return genome.chromosomes[chr][index].end(), genome.chromosomes[chr][index+1].start(),'+'
                else:
                    return None,None,None
            elif orient1 == 't' and index > 0:
                if genome.chromosomes[chr][index-1].name() == m2:
                    if orient2 == 't' and genome.chromosomes[chr][index-1].orientation() < 0:
                        return genome.chromosomes[chr][index-1].end(), genome.chromosomes[chr][index].start(),'-'
                    elif orient2 == 'h' and genome.chromosomes[chr][index-1].orientation() > 0:
                        return genome.chromosomes[chr][index-1].end(), genome.chromosomes[chr][index].start(),'-'
                else:
                    return None,None,None
        elif orientation < 0:
            if orient1 == 'h' and index > 0:
                if genome.chromosomes[chr][index-1].name() == m2:
                    if orient2 == 't' and genome.chromosomes[chr][index-1].orientation() < 0:
                        return genome.chromosomes[chr][index-1].end(), genome.chromosomes[chr][index].start(),'-'
                    elif orient2 == 'h' and genome.chromosomes[chr][index-1].orientation() > 0 :
                        return genome.chromosomes[chr][index-1].end(), genome.chromosomes[chr][index].start(),'-'
                else:
                    return None,None,None
            elif orient1 == 't' and index < len(genome.chromosomes[chr])-1:
                if genome.chromosomes[chr][index+1].name() == m2:
                    if orient2 == 't' and genome.chromosomes[chr][index+1].orientation() > 0:
                        return genome.chromosomes[chr][index].end(), genome.chromosomes[chr][index+1].start(),'+'
                    elif orient2 == 'h' and genome.chromosomes[chr][index+1].orientation() < 0:
                        return genome.chromosomes[chr][index].end(), genome.chromosomes[chr][index+1].start(),'+'
                else:
                    return None,None,None

    return None,None,None

def doubler(string):
    num = abs(int(string))
    tail = str(2*num-1)
    head = str(2*num)
    if string[0] == '-':
        return head+' '+tail
    else:
        return tail+' '+head

o = open(sys.argv[5],'w')
for adj in range(len(c1p_adjacencies)):
    #print c1p_adjacencies[adj]
#    adj1,adj2 = c1p_adjacencies[adj]
    adj1,adj2 = sorted(c1p_adjacencies[adj])
    if adj1 == makeMate(adj2):
        continue
    # Rubberducking:
    # If adj1,adj2 == (2,3),(3,2), then both m1,m2 are positively oriented, i.e. lexicographically.
    # If adj1,adj2 == (2,4),(4,2), then m1 is +, m2 is -.
    # If adj1,adj2 == (1,3),(3,1), then m1 is -, m2 is +.
    # If adj1,adj2 == (1,4),(4,1), then m1 is -, m2 is -.
    # Relabelled using _h,_t to make it clear, since markers are not doubled.
    o.write(">"+adj1+' '+adj2+"\n")
    m1 = int(adj1)
    m2 = int(adj2)
    orient1 = 'h'
    orient2 = 't'
    if m1%2 == 1:
        m1 += 1
        orient1 = 't'
    if m2%2 == 0:
        orient2 = 'h'
    else:
        m2 += 1
    m1 = str(m1>>1)
    m2 = str(m2>>1)
    for loc1 in synteny_blocks_big[m1].loci:
        if loc1.species not in all_species or loc1.species not in sp_list[adj]:
            continue
        gap_data = findGapPositions(genomes_big[loc1.species],m1,orient1,m2,orient2,loc1)
        if gap_data[0] != None:
            a = (adj1,adj2)
            adj_position = (loc1.species,loc1.chromosome,gap_data[0],gap_data[1],gap_data[2])
            temp = getContainment(genomes_small[loc1.species],loc1.chromosome,gap_data[0],gap_data[1])
            o.write(loc1.species+'.'+loc1.chromosome+':')
            new_temp = [doubler(x) for x in temp]
            o.write(' '.join(new_temp)+' ')
            # o.write(big_orientation+"\n")
            o.write(gap_data[-1]+"\n")
        else:
            continue
    o.write("\n")
    if adj == intercar_pos:
        o.write("#INTERCAR\n")

o.close()
#o = open("../Mammalian_data/MGRA_comparison/blocks_adjacencies_20_in_100",'w')
#for adj in range(1,len(c1p_adjacencies)+1):
#    adj1,adj2 = c1p_adjacencies[adj]
#    # Rubberducking:
#    # If adj1,adj2 == (2,3),(3,2), then both m1,m2 are positively oriented, i.e. lexicographically.
#    # If adj1,adj2 == (2,4),(4,2), then m1 is +, m2 is -.
#    # If adj1,adj2 == (1,3),(3,1), then m1 is -, m2 is +.
#    # If adj1,adj2 == (1,4),(4,1), then m1 is -, m2 is -.
#    # Relabelled using _h,_t to make it clear, since markers are not doubled.
#    o.write(">"+adj1[0]+'_'+adj1[1]+' '+adj2[0]+'_'+adj2[1]+"\n")
#    m1,orient1 = adj1
#    m2,orient2 = adj2
#    for loc1 in uu_synteny_blocks_100k[m1].loci:
#        gap_data = findGapPositions(uu_genomes_100k[loc1.species],m1,orient1,m2,orient2,loc1)
#        if gap_data[0] != None:
#            a = (adj1,adj2)
#            adj_position = (loc1.species,loc1.chromosome,gap_data[0],gap_data[1],gap_data[2])
#            temp = getContainment(genomes_20k[loc1.species],loc1.chromosome,gap_data[0],gap_data[1])
#            o.write(loc1.species+'.'+loc1.chromosome+':')
#            for t in temp:
#                o.write(t+' ')
#            o.write(gap_data[-1]+"\n")
#        else:
#            continue
#    o.write("\n")
#    if index == intercar_pos:
#        o.write("#INTERCAR\n")
#
#o.close()
#o = open("../Mammalian_data/MGRA_comparison/blocks_adjacencies_100_in_500",'w')
#for index in c1p_adjacencies:
#    adj1,adj2 = c1p_adjacencies[index]
#    # Rubberducking:
#    # If adj1,adj2 == (2,3),(3,2), then both m1,m2 are positively oriented, i.e. lexicographically.
#    # If adj1,adj2 == (2,4),(4,2), then m1 is +, m2 is -.
#    # If adj1,adj2 == (1,3),(3,1), then m1 is -, m2 is +.
#    # If adj1,adj2 == (1,4),(4,1), then m1 is -, m2 is -.
#    # Relabelled using _h,_t to make it clear, since markers are not doubled.
#    o.write(">"+adj1[0]+'_'+adj1[1]+' '+adj2[0]+'_'+adj2[1]+"\n")
#    m1,orient1 = adj1
#    m2,orient2 = adj2
#    for loc1 in uu_synteny_blocks_500k[m1].loci:
#        gap_data = findGapPositions(uu_genomes_500k[loc1.species],m1,orient1,m2,orient2,loc1)
#        if gap_data[0] != None:
#            a = (adj1,adj2)
#            adj_position = (loc1.species,loc1.chromosome,gap_data[0],gap_data[1],gap_data[2])
#            temp = getContainment(genomes_100k[loc1.species],loc1.chromosome,gap_data[0],gap_data[1])
#            o.write(loc1.species+'.'+loc1.chromosome+':')
#            for t in temp:
#                o.write(t+' ')
#            o.write(gap_data[-1]+"\n")
#        else:
#            continue
#    o.write("\n")
#    if index == intercar_pos:
#        o.write("#INTERCAR\n")
#o.close()
#o = open("../Mammalian_data/MGRA_comparison/blocks_adjacencies_300_in_500",'w')
#for index in c1p_adjacencies:
#    adj1,adj2 = c1p_adjacencies[index]
#    # Rubberducking:
#    # If adj1,adj2 == (2,3),(3,2), then both m1,m2 are positively oriented, i.e. lexicographically.
#    # If adj1,adj2 == (2,4),(4,2), then m1 is +, m2 is -.
#    # If adj1,adj2 == (1,3),(3,1), then m1 is -, m2 is +.
#    # If adj1,adj2 == (1,4),(4,1), then m1 is -, m2 is -.
#    # Relabelled using _h,_t to make it clear, since markers are not doubled.
#    o.write(">"+adj1[0]+'_'+adj1[1]+' '+adj2[0]+'_'+adj2[1]+"\n")
#    m1,orient1 = adj1
#    m2,orient2 = adj2
#    for loc1 in uu_synteny_blocks_500k[m1].loci:
#        gap_data = findGapPositions(uu_genomes_500k[loc1.species],m1,orient1,m2,orient2,loc1)
#        if gap_data[0] != None:
#            a = (adj1,adj2)
#            adj_position = (loc1.species,loc1.chromosome,gap_data[0],gap_data[1],gap_data[2])
#            temp = getContainment(genomes_300k[loc1.species],loc1.chromosome,gap_data[0],gap_data[1])
#            o.write(loc1.species+'.'+loc1.chromosome+':')
#            for t in temp:
#                o.write(t+' ')
#            o.write(gap_data[-1]+"\n")
#        else:
#            continue
#    o.write("\n")
#    if index == intercar_pos:
#        o.write("#INTERCAR\n")
#
#o.close()
#o = open("../Mammalian_data/MGRA_comparison/blocks_adjacencies_50_in_300",'w')
#for index in c1p_adjacencies:
#    adj1,adj2 = c1p_adjacencies[index]
#    # Rubberducking:
#    # If adj1,adj2 == (2,3),(3,2), then both m1,m2 are positively oriented, i.e. lexicographically.
#    # If adj1,adj2 == (2,4),(4,2), then m1 is +, m2 is -.
#    # If adj1,adj2 == (1,3),(3,1), then m1 is -, m2 is +.
#    # If adj1,adj2 == (1,4),(4,1), then m1 is -, m2 is -.
#    # Relabelled using _h,_t to make it clear, since markers are not doubled.
#    o.write(">"+adj1[0]+'_'+adj1[1]+' '+adj2[0]+'_'+adj2[1]+"\n")
#    m1,orient1 = adj1
#    m2,orient2 = adj2
#    for loc1 in uu_synteny_blocks_300k[m1].loci:
#        gap_data = findGapPositions(uu_genomes_300k[loc1.species],m1,orient1,m2,orient2,loc1)
#        if gap_data[0] != None:
#            a = (adj1,adj2)
#            adj_position = (loc1.species,loc1.chromosome,gap_data[0],gap_data[1],gap_data[2])
#            temp = getContainment(genomes_50k[loc1.species],loc1.chromosome,gap_data[0],gap_data[1])
#            o.write(loc1.species+'.'+loc1.chromosome+':')
#            for t in temp:
#                o.write(t+' ')
#            o.write(gap_data[-1]+"\n")
#        else:
#            continue
#    o.write("\n")
#    if index == intercar_pos:
#        o.write("#INTERCAR\n")
#
#o.close()
#o = open("../Mammalian_data/MGRA_comparison/blocks_adjacencies_20_in_300",'w')
#for index in c1p_adjacencies:
#    adj1,adj2 = c1p_adjacencies[index]
#    # Rubberducking:
#    # If adj1,adj2 == (2,3),(3,2), then both m1,m2 are positively oriented, i.e. lexicographically.
#    # If adj1,adj2 == (2,4),(4,2), then m1 is +, m2 is -.
#    # If adj1,adj2 == (1,3),(3,1), then m1 is -, m2 is +.
#    # If adj1,adj2 == (1,4),(4,1), then m1 is -, m2 is -.
#    # Relabelled using _h,_t to make it clear, since markers are not doubled.
#    o.write(">"+adj1[0]+'_'+adj1[1]+' '+adj2[0]+'_'+adj2[1]+"\n")
#    m1,orient1 = adj1
#    m2,orient2 = adj2
#    for loc1 in uu_synteny_blocks_300k[m1].loci:
#        gap_data = findGapPositions(uu_genomes_300k[loc1.species],m1,orient1,m2,orient2,loc1)
#        if gap_data[0] != None:
#            a = (adj1,adj2)
#            adj_position = (loc1.species,loc1.chromosome,gap_data[0],gap_data[1],gap_data[2])
#            temp = getContainment(genomes_20k[loc1.species],loc1.chromosome,gap_data[0],gap_data[1])
#            o.write(loc1.species+'.'+loc1.chromosome+':')
#            for t in temp:
#                o.write(t+' ')
#            o.write(gap_data[-1]+"\n")
#        else:
#            continue
#    o.write("\n")
#    if index == intercar_pos:
#        o.write("#INTERCAR\n")
#
#o.close()
#o = open("../Mammalian_data/MGRA_comparison/blocks_adjacencies_50_in_500",'w')
#for index in c1p_adjacencies:
#    adj1,adj2 = c1p_adjacencies[index]
#    # Rubberducking:
#    # If adj1,adj2 == (2,3),(3,2), then both m1,m2 are positively oriented, i.e. lexicographically.
#    # If adj1,adj2 == (2,4),(4,2), then m1 is +, m2 is -.
#    # If adj1,adj2 == (1,3),(3,1), then m1 is -, m2 is +.
#    # If adj1,adj2 == (1,4),(4,1), then m1 is -, m2 is -.
#    # Relabelled using _h,_t to make it clear, since markers are not doubled.
#    o.write(">"+adj1[0]+'_'+adj1[1]+' '+adj2[0]+'_'+adj2[1]+"\n")
#    m1,orient1 = adj1
#    m2,orient2 = adj2
#    for loc1 in uu_synteny_blocks_500k[m1].loci:
#        gap_data = findGapPositions(uu_genomes_500k[loc1.species],m1,orient1,m2,orient2,loc1)
#        if gap_data[0] != None:
#            a = (adj1,adj2)
#            adj_position = (loc1.species,loc1.chromosome,gap_data[0],gap_data[1],gap_data[2])
#            temp = getContainment(genomes_50k[loc1.species],loc1.chromosome,gap_data[0],gap_data[1])
#            o.write(loc1.species+'.'+loc1.chromosome+':')
#            for t in temp:
#                o.write(t+' ')
#            o.write(gap_data[-1]+"\n")
#        else:
#            continue
#    o.write("\n")
#    if index == intercar_pos:
#        o.write("#INTERCAR\n")
#
#o.close()
#o = open("../Mammalian_data/MGRA_comparison/blocks_adjacencies_20_in_500",'w')
#for index in c1p_adjacencies:
#    adj1,adj2 = c1p_adjacencies[index]
#    # Rubberducking:
#    # If adj1,adj2 == (2,3),(3,2), then both m1,m2 are positively oriented, i.e. lexicographically.
#    # If adj1,adj2 == (2,4),(4,2), then m1 is +, m2 is -.
#    # If adj1,adj2 == (1,3),(3,1), then m1 is -, m2 is +.
#    # If adj1,adj2 == (1,4),(4,1), then m1 is -, m2 is -.
#    # Relabelled using _h,_t to make it clear, since markers are not doubled.
#    o.write(">"+adj1[0]+'_'+adj1[1]+' '+adj2[0]+'_'+adj2[1]+"\n")
#    m1,orient1 = adj1
#    m2,orient2 = adj2
#    for loc1 in uu_synteny_blocks_500k[m1].loci:
#        gap_data = findGapPositions(uu_genomes_500k[loc1.species],m1,orient1,m2,orient2,loc1)
#        if gap_data[0] != None:
#            a = (adj1,adj2)
#            adj_position = (loc1.species,loc1.chromosome,gap_data[0],gap_data[1],gap_data[2])
#            temp = getContainment(genomes_20k[loc1.species],loc1.chromosome,gap_data[0],gap_data[1])
#            o.write(loc1.species+'.'+loc1.chromosome+':')
#            for t in temp:
#                o.write(t+' ')
#            o.write(gap_data[-1]+"\n")
#        else:
#            continue
#    o.write("\n")
#    if index == intercar_pos:
#        o.write("#INTERCAR\n")
#o.close()
#
