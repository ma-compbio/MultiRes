#!/usr/bin/env

import sys
import markers
import bisect

# Genome data structure, kept ordered, allows in-order insertions
class Genome:
    # Constructor
    def __init__(self,sp):
        self.species = sp
        self.chromosomes = {}
    # Add locus at appropriate location to correct chromosome
    def addLocus(self,mark,locus,threshold=0):
        if not(abs(locus.start-locus.end)+1 < threshold):
            chrom = locus.chromosome
            if not chrom in self.chromosomes:
                self.chromosomes[chrom] = []
            i = bisect.bisect_left([x.start() for x in self.chromosomes[chrom]],locus.start)
            self.chromosomes[chrom].insert(i,Marker(mark,locus.start,locus.end,locus.orientation))
    # Find index of a marker in a chromosome, given its start location in the chromosome.
    def findLocus(self,chrom,location):
        current_chrom = self.chromosomes[chrom]
        i = bisect.bisect_left([x.start() for x in current_chrom],location)
        if current_chrom[i].start() == location:
            return i
        else:
            return -sys.maxint+1 
    def length(self,chrom):
        return len(self.chromosomes[chrom])
    # Given a location, either finds marker starting at that location, or first marker overlapping with it.
    def findLocusFuzzy(self,chrom,location):
        try:
            current_chrom = self.chromosomes[chrom]
            i = bisect.bisect_left([x.start() for x in current_chrom],location)
            if i < len(current_chrom):
                if current_chrom[i-1].end() >= location and current_chrom[i-1].start() <= location:
                    return i-1
                elif current_chrom[i].start() >= location:
                    return i
            elif i == len(current_chrom):
                return -sys.maxint+1
        except:
            return -sys.maxint+1
    # Reverse a segment of the genome        
    def reverseSubChrom(self,chrom,start=0,end=-1):
        if end < 0:
            end = len(self.chromosomes[chrom])
        temporary_segment = self.chromosomes[chrom][start:end][:]
        start_pos = self.chromosomes[chrom][start].start()
        temporary_segment = temporary_segment[::-1]
        k = 0
        for i in range(start,end):
            temp_len = temporary_segment[k].end()-temporary_segment[k].start()
            self.chromosomes[chrom][i] = Marker(temporary_segment[k].name(),\
                    start_pos,\
                    start_pos+temp_len,\
                    -1*temporary_segment[k].orientation())
            if k < len(temporary_segment)-1:
                prev_pos = start_pos
                start_pos += temp_len - temporary_segment[k+1].end()+\
                        temporary_segment[k].start()
            k += 1
            
    # Cut a segment and insert at another position            
    def translocate(self,chrom0,chrom1,pos00,pos01,pos1,orient):
        if chrom0 == chrom1:
            if pos1 < pos00:
                temp_list = self.chromosomes[chrom0][pos00:pos01][:]
                self.deleteFromChrom(chrom0,pos00,pos01)
                self.insertToChrom(chrom1,pos1,temp_list,orient)
                del temp_list
            elif pos1 > pos01:
                temp_list = self.chromosomes[chrom0][pos00:pos01][:]
                self.insertToChrom(chrom1,pos1,temp_list,orient)
                self.deleteFromChrom(chrom0,pos00,pos01)
                del temp_list
        else:
            temp_list = self.chromosomes[chrom0][pos00:pos01][:]
            self.insertToChrom(chrom1,pos1,temp_list,orient)
            self.deleteFromChrom(chrom0,pos00,pos01)
            del temp_list
    # Insert a segment in a given orientation within a given chromosome            
    def insertToChrom(self,chrom,pos,marker_list,orient):
        segment_len = marker_list[-1].end()-marker_list[0].start()+1
        if pos > 0:
            if pos < len(self.chromosomes[chrom]):
                gap_len = (self.chromosomes[chrom][pos].start()-\
                        self.chromosomes[chrom][pos-1].end())>>1
                start_pos = self.chromosomes[chrom][pos-1].end()+gap_len
            else:
                gap_len = 0
                start_pos = self.chromosomes[chrom][-1].end()+1
        else:
            gap_len = 0
            start_pos = 1
        for i in range(pos,len(self.chromosomes[chrom])):
            self.chromosomes[chrom][i]._start += segment_len - gap_len
            self.chromosomes[chrom][i]._end += segment_len - gap_len
        if orient > 0:
            current_pos = pos
            for i in range(len(marker_list)):
                m = marker_list[i]
                m_len = m.end()-m.start()+1
                self.chromosomes[chrom].insert(current_pos,\
                        Marker(m.name(),start_pos,start_pos+m_len,m.orientation()))
                if i < len(marker_list)-1:
                    start_pos += m_len + marker_list[i+1].start()-m.end()
                current_pos += 1                    
        elif orient < 0:
            current_pos = pos
            for i in range(len(marker_list)-1,-1,-1):
                m = marker_list[i]
                m_len = abs(m.end()-m.start())+1
                self.chromosomes[chrom].insert(current_pos,\
                        Marker(m.name(),start_pos,start_pos+m_len,-1*m.orientation()))
                if i > 0:
                    start_pos += m_len + m.start()-marker_list[i-1].end()
                current_pos += 1
    # Delete range of markers from a chromosome                
    def deleteFromChrom(self,chrom,start,end):
        seg_len = self.chromosomes[chrom][end-1].end()-self.chromosomes[chrom][start].start()
        for i in range(end,len(self.chromosomes[chrom])):
            self.chromosomes[chrom][i]._start -= seg_len
            self.chromosomes[chrom][i]._end -= seg_len
        self.chromosomes[chrom] = self.chromosomes[chrom][:start]+self.chromosomes[chrom][end:]
        
def resizeGenome(genome,chrom,slice_list,start,end=-1):
    start_pos = genome.chromosomes[chrom][start][0]
    for i in range(len(slice_list)):
        l = slice_list[i][1]-slice_list[i][0]
        self.positions[chrom][start] = (start_pos,start_pos+l)
        if i < len(slice_list)-1:
            gap = slice_list[i+1][0] - slice_list[i][1]
            start_pos = self.positions[chrom][start][1]
        start += 1
        
# Simplified data structure with marker id + positions + orientation
class Marker:
    def __init__(self,name,start,end,orient):
        self._name = name
        self._start = start
        self._end = end
        self._orientation = orient
    def name(self):
        return self._name
    def start(self):
        return self._start
    def end(self):
        return self._end
    def length(self):
        return abs(self.start()-self.end())+1
    def orientation(self):
        return self._orientation


# Get dictionary of sorted genomes
def getGenomes(all_markers,sp,threshold=-1):
    # Dictionary of genomes
    genomes = {}
    # Create for each species in list
    for s in sp:
        genomes[s] = Genome(s)
    # Add markers in place
    for fam in all_markers:
        try:
            for locus in all_markers[fam].loci:
                if locus.species in sp:
                    genomes[locus.species].addLocus(all_markers[fam].id,locus,threshold)
        except TypeError:
            for locus in fam.loci:
                if locus.species in sp:
                    genomes[locus.species].addLocus(fam.id,locus,threshold)

    return genomes
