#!/usr/bin/env
import sys
from DataStructures import genomes
from DataStructures import markers

species_list_in = []
species_list_out = []
with open(sys.argv[3]) as sp_list_file:
    toggle_type = False
    for line in sp_list_file:
        if "INGROUPS" in line.upper():
            toggle_type = True
        elif "OUTGROUPS" in line.upper():
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
#synteny_blocks_20k = markers.read_hom_families_file("../Mammalian_data/Conserved.Segments.20K.Species_filtered.uu")
#synteny_blocks_50k = markers.read_hom_families_file("../Mammalian_data/Conserved.Segments.50K.Species_filtered")

all_species = species_list_in+species_list_out

#uu_genomes_100k = genomes.getGenomes(uu_synteny_blocks_100k,all_species)
#uu_genomes_300k = genomes.getGenomes(uu_synteny_blocks_300k,all_species)
#uu_genomes_500k = genomes.getGenomes(uu_synteny_blocks_500k,all_species)
#uu_genomes_20k = genomes.getGenomes(uu_synteny_blocks_20k,all_species)
#uu_genomes_50k = genomes.getGenomes(uu_synteny_blocks_50k,all_species)

#genomes_100k = genomes.getGenomes(synteny_blocks_100k,all_species)
#genomes_300k = genomes.getGenomes(synteny_blocks_300k,all_species)
#genomes_500k = genomes.getGenomes(synteny_blocks_500k,all_species)
#genomes_20k = genomes.getGenomes(synteny_blocks_20k,all_species)
#genomes_50k = genomes.getGenomes(synteny_blocks_50k,all_species)
#genomes_big = genomes.getGenomes(synteny_blocks_big,all_species)
genomes_small  = genomes.getGenomes(synteny_blocks_small,all_species)

def getContainment(genome,chromosome,start,end):
    index = genome.findLocusFuzzy(chromosome,start)
    # print genome.chromosomes[chromosome][index].name(),genome.chromosomes[chromosome][index].start(),genome.chromosomes[chromosome][index].end()
    temp_list = []
    first_marker = True
    while index != -sys.maxint+1 and index < len(genome.chromosomes[chromosome]) and genome.chromosomes[chromosome][index].start() < end:
        orient = ''
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

def doubler(string):
    num = abs(int(string))
    tail = str(2*num-1)
    head = str(2*num)
    if string[0] == '-':
        return head+' '+tail
    else:
        return tail+' '+head


o = open(sys.argv[4],'w')
for mbig in synteny_blocks_big:
    o.write(">"+mbig+"\n")
    for locus in synteny_blocks_big[mbig].loci:
        sp = locus.species
        if sp not in all_species:
            continue
        chrom = locus.chromosome
        big_start = locus.start
        big_end = locus.end
        big_orientation = ''
        if locus.orientation < 0:
            big_orientation = '-'
        elif locus.orientation > 0:
            big_orientation = '+'
        o.write(sp+"."+chrom+": ")
        temp = getContainment(genomes_small[sp],chrom,big_start,big_end)
#        new_temp = temp#[doubler(x) for x in temp]
        new_temp = [doubler(x) for x in temp]
        o.write(' '.join(new_temp)+' ')
        o.write(big_orientation+"\n")
    o.write("\n")
o.close()

#o = open("../Mammalian_data/blocks_100_in_500_sp_filtered",'w')
#for m500 in uu_synteny_blocks_500k:
#    o.write(">"+m500+"\n")
#    for locus in uu_synteny_blocks_500k[m500].loci:
#        sp = locus.species
#        chrom = locus.chromosome
#        big_start = locus.start
#        big_end = locus.end
#        big_orientation = ''
#        if locus.orientation < 0:
#            big_orientation = '-'
#        elif locus.orientation > 0:
#            big_orientation = '+'
#        o.write(sp+"."+chrom+": ")
#        temp = getContainment(genomes_100k[sp],chrom,big_start,big_end)
#        for t in temp:
#            o.write(t+" ")
#        o.write(big_orientation+"\n")
#    o.write("\n")
#o.close()

#o = open("../Mammalian_data/MGRA_comparison/blocks_300_in_500_sp_filtered",'w')
#for m500 in uu_synteny_blocks_500k:
#    o.write(">"+m500+"\n")
#    for locus in uu_synteny_blocks_500k[m500].loci:
#        sp = locus.species
#        chrom = locus.chromosome
#        big_start = locus.start
#        big_end = locus.end
#        big_orientation = ''
#        if locus.orientation < 0:
#            big_orientation = '-'
#        elif locus.orientation > 0:
#            big_orientation = '+'
#        o.write(sp+"."+chrom+": ")
#        temp = getContainment(genomes_300k[sp],chrom,big_start,big_end)
#        for t in temp:
#            o.write(t+" ")
#        o.write(big_orientation+"\n")
#    o.write("\n")
#o.close()

#o = open("../Mammalian_data/blocks_100_in_300_sp_filtered",'w')
#for m300 in uu_synteny_blocks_300k:
#    o.write(">"+m300+"\n")
#    for locus in uu_synteny_blocks_300k[m300].loci:
#        sp = locus.species
#        chrom = locus.chromosome
#        big_start = locus.start
#        big_end = locus.end
#        big_orientation = ''
#        if locus.orientation < 0:
#            big_orientation = '-'
#        elif locus.orientation > 0:
#            big_orientation = '+'
#        o.write(sp+"."+chrom+": ")
#        temp = getContainment(genomes_100k[sp],chrom,big_start,big_end)
#        for t in temp:
#            o.write(t+" ")
#        o.write(big_orientation+"\n")
#    o.write("\n")
#o.close()
#

#o = open("../Mammalian_data/MGRA_comparison/blocks_50_in_500_sp_filtered",'w')
#for m500 in uu_synteny_blocks_500k:
#    o.write(">"+m500+"\n")
#    for locus in uu_synteny_blocks_500k[m500].loci:
#        sp = locus.species
#        chrom = locus.chromosome
#        big_start = locus.start
#        big_end = locus.end
#        big_orientation = ''
#        if locus.orientation < 0:
#            big_orientation = '-'
#        elif locus.orientation > 0:
#            big_orientation = '+'
#        o.write(sp+"."+chrom+": ")
#        temp = getContainment(genomes_50k[sp],chrom,big_start,big_end)
#        for t in temp:
#            o.write(t+" ")
#        o.write(big_orientation+"\n")
#    o.write("\n")
#o.close()
#
#o = open("../Mammalian_data/MGRA_comparison/blocks_20_in_500_sp_filtered",'w')
#for m500 in uu_synteny_blocks_500k:
#    o.write(">"+m500+"\n")
#    for locus in uu_synteny_blocks_500k[m500].loci:
#        sp = locus.species
#        chrom = locus.chromosome
#        big_start = locus.start
#        big_end = locus.end
#        big_orientation = ''
#        if locus.orientation < 0:
#            big_orientation = '-'
#        elif locus.orientation > 0:
#            big_orientation = '+'
#        o.write(sp+"."+chrom+": ")
#        temp = getContainment(genomes_20k[sp],chrom,big_start,big_end)
#        for t in temp:
#            o.write(t+" ")
#        o.write(big_orientation+"\n")
#    o.write("\n")
#o.close()
#

#o = open("../Mammalian_data/MGRA_comparison/blocks_20_in_300_sp_filtered",'w')
#for m300 in uu_synteny_blocks_300k:
#    o.write(">"+m300+"\n")
#    for locus in uu_synteny_blocks_300k[m300].loci:
#        sp = locus.species
#        chrom = locus.chromosome
#        big_start = locus.start
#        big_end = locus.end
#        big_orientation = ''
#        if locus.orientation < 0:
#            big_orientation = '-'
#        elif locus.orientation > 0:
#            big_orientation = '+'
#        o.write(sp+"."+chrom+": ")
#        temp = getContainment(genomes_20k[sp],chrom,big_start,big_end)
#        for t in temp:
#            o.write(t+" ")
#        o.write(big_orientation+"\n")
#    o.write("\n")
#o.close()
#
#o = open("../Mammalian_data/MGRA_comparison/blocks_50_in_300_sp_filtered",'w')
#for m300 in uu_synteny_blocks_300k:
#    o.write(">"+m300+"\n")
#    for locus in uu_synteny_blocks_300k[m300].loci:
#        sp = locus.species
#        chrom = locus.chromosome
#        big_start = locus.start
#        big_end = locus.end
#        big_orientation = ''
#        if locus.orientation < 0:
#            big_orientation = '-'
#        elif locus.orientation > 0:
#            big_orientation = '+'
#        o.write(sp+"."+chrom+": ")
#        temp = getContainment(genomes_50k[sp],chrom,big_start,big_end)
#        for t in temp:
#            o.write(t+" ")
#        o.write(big_orientation+"\n")
#    o.write("\n")
#o.close()
#
#o = open("../Mammalian_data/MGRA_comparison/blocks_20_in_100_sp_filtered",'w')
#for m100 in uu_synteny_blocks_100k:
#    o.write(">"+m100+"\n")
#    for locus in uu_synteny_blocks_100k[m100].loci:
#        sp = locus.species
#        chrom = locus.chromosome
#        big_start = locus.start
#        big_end = locus.end
#        big_orientation = ''
#        if locus.orientation < 0:
#            big_orientation = '-'
#        elif locus.orientation > 0:
#            big_orientation = '+'
#        o.write(sp+"."+chrom+": ")
#        temp = getContainment(genomes_20k[sp],chrom,big_start,big_end)
#        for t in temp:
#            o.write(t+" ")
#        o.write(big_orientation+"\n")
#    o.write("\n")
#o.close()
#
