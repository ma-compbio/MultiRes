#!/usr/bin/env
import sys

# locus
class Locus:
    # Locus constructor.
    # Parameters
    # species: str - species name
    # chromosome: str - chromosome name
    # start: int - locus start
    # end: int - locus end
    # orientation: int - locus orientation (positive for forwards, negative for reversed and 0 for unoriented)
    # comment: str - comment
    def __init__(self, species, chromosome, start, end, orientation, comment=''):
        self.species = species
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.orientation = orientation
        self.comment = comment
    def length(self):
        return abs(self.end-self.start)+1
    def doubler(self):
        end0 = Locus(self.species,self.chromosome,self.start,self.start,0)
        end1 = Locus(self.species,self.chromosome,self.end,self.end,0)
        if self.orientation > 0:
            return end0,end1,True
        elif self.orientation < 0:
            return end1,end0,True
        else:
            return end0,end1,False
    # Compares two loci. Compares species name then chromosome name then locus start.
    # Parameters
    # other: Locus
    # Return: int - compare value
    def __cmp__(self, other):
        if type(other) is Locus:
            a = cmp(self.species, other.species)
            b = cmp(self.chromosome, other.chromosome)      
            if a == 0:
                if b == 0:
                    return cmp(self.start, other.start)
                return b
            return a
        else:
            return 1

    # Returns a string representing the locus with given the separators. Will not print the comment is it is empty.
    # <species name>sp1<chromosome name>sp2<start>sp3<end>sp4<orientation>sp5<comment>
    # Parameters
    # sp1: str
    # sp2: str
    # sp3: str
    # sp4: str
    # sp5: str
    # Return: str - string representation of the locus
    def str_general(self,sp1,sp2,sp3,sp4,sp5):
        if self.orientation > 0:
            orient_str = "+"
        elif self.orientation < 0:
            orient_str = "-"
        else: # == 0
            orient_str = "X"
        #endif
        if len(self.comment) > 0:
            comment_str = sp5 +  self.comment
        else:
            comment_str = ""
        #endif
        return self.species + sp1 + self.chromosome + sp2 + str(self.start) + sp3 + str(self.end) + sp4 +orient_str + comment_str
    #enddef

    # Returns a string representing the locus. Will not print the comment is it is empty.
    # <species name>.<chromosome name>:<start>-<end> <orientation> #<comment>
    # or (if the species name contains '.')
    # <species name> <chromosome name> <start> <end> <orientation> #<comment>
    # Return: str - string representation of the locus
    def __str__(self):
        if self.species.find(".") < 0:
            return self.str_general(".",":","-"," "," #")
        else:
            return self.str_general(" "," "," "," "," #")
        #endif
    #enddef
    
    # Static method which parses a locus from a string.
    # <species name>.<chromosome name>:<start>-<end> <orientation> #<comment>
    # or
    # <species name> <chromosome name> <start> <end> <orientation> #<comment>
    # Parameters
    # string - str: the string to parse
    # Return Locus - the locus in the file (None on failure)
    @staticmethod
    def from_string(string,threshold = -1):
        #assert type(string) is str
    
        trunc_str = string.strip()
        full_str = trunc_str
    
        if len(trunc_str) == 0 or trunc_str[0] == '>' or trunc_str[0] == '#':
            return None
        #endif

        # Extracting the comments field
        comment_split = trunc_str.split('#', 1)
        if len(comment_split) > 1:
            comment = comment_split[1]
        else:
            comment = ""
        #endif
        trunc_str = comment_split[0].strip()
        
        # using space separated format
        split_str =  trunc_str.split()
        if len(split_str) >= 4:
            species = split_str[0]
            chromosome = split_str[1]
            
            try:
                start = int(split_str[2])       
            except:
                print("Warning: start coordinate not an integer for locus. String : \'" + full_str + "\'")
                return None
            #endtry
            
            try:
                end = int(split_str[2])     
            except:
                print("Warning: end coordinate not an integer for locus. String : \'" + full_str + "\'")
                return None
            #endtry
            
            if len(split_str) == 4:
                orientation = 0
            else:
                trunc_str = split_str[4].strip()
                print trunc_str
                if trunc_str[0] == '+':
                    orientation = 1         
                elif trunc_str[0] == '-':
                    orientation = -1
                elif trunc_str[0].lower() == 'x':
                    orientation = 0
                else:
                    print("Warning: unknown orientation for locus. String : \'" + full_str + "\'")
                    return None
                #endif
            #endif
            if abs(start-end)+1 < threshold:
                return None
            return Locus(species, chromosome, start, end, orientation, comment)
        #endif
        
        # using old format      
        species, trunc_str = Locus.obj_split(trunc_str, ".", 2, full_str, "species")
        if trunc_str == None:
            return None
        #endif
        
        chromosome, trunc_str = Locus.obj_split(trunc_str, ":", 2, full_str, "chromosome")
        if trunc_str == None:
            return None
        #endif
        
        start, trunc_str = Locus.obj_split(trunc_str, "-", 2, full_str, "start")
        if trunc_str == None:
            return None
        #endif
        try:
            start = int(start)      
        except:
            print("Warning: start coordinate not an integer for locus. String : \'" + full_str + "\'")
            return None
        #endtry
                
        end, trunc_str = Locus.obj_split(trunc_str, " ", 1, full_str, "end")
        if trunc_str == None:
            return None
        #endif
        try:
            end = int(end)      
        except:
            print("Warning: end coordinate not an integer for locus. String : \'" + full_str + "\'")
            return None
        #endtry
        
        if len(trunc_str) == 0:
            orientation = 0
        else:
            trunc_str = trunc_str.strip()
            if trunc_str[0] == '+':
                orientation = 1         
            elif trunc_str[0] == '-':
                orientation = -1
            elif trunc_str[0].lower() == 'x':
                orientation = 0
            else:
                print("Warning: unknown orientation for locus. String : \'" + full_str + "\'")
                return None
            #endif
        #endif
        if abs(start-end)+1 < threshold:
            return None
        return Locus(species, chromosome, start, end, orientation, comment)
    #enddef
    
    # Static auxilary method for from_string.
    # Parameters:
    # string - str
    # delimeter - str
    # size - int
    # full_string - str
    # warning name - str
    # Return - str, str
    @staticmethod
    def obj_split(string, delimiter, size, full_string, warning_name):
        #assert type(string) is str
        #assert type(delimiter) is str
        #assert type(size) is int
        #assert type(full_string) is str
        #assert type(warning_name) is str
    
        split = string.split(delimiter, 1)
                
        if len(split) < size:
            print("Warning: could not process " + warning_name + " for locus. String: \'" + full_string + "\'")
            
            return None, None
        #endif
                    
        if len(split) < 2:
            split.append(None)
        #endif
        
        return split[0], split[1]
    #enddef 
#endclass

# homologous family
class HomFam:
    # HomFam constructor
    # ident - str: marker identity
    # loci - list of Locus: hom. family loci
    # copy_number - int: copy number
    # comment - str: comment
    def __init__(self, ident, loci, copy_number, comment=''):
        #assert type(ident) is str
        #assert type(loci) is list
        #assert type(copy_number) is int
        #assert type(comment) is str
    
        self.id = ident
        self.loci = loci
        self.copy_number = copy_number
        self.comment = comment
    #enddef
    def doubler(self,num=0):
        if num == 1:
            head = str(2*int(self.id)-1)
            tail = str(2*int(self.id))
        else:            
            head = self.id+'_h'
            tail = self.id+'_t'
        copy_n = self.copy_number
        new_comment = self.comment
        head_loci = []
        tail_loci = []
        for l in self.loci:
            head_locus,tail_locus, doubleable = l.doubler()
            if not(doubleable):
                print "Error: Markers have no orientation."
                sys.exit(-1)
            else:
                head_loci.append(head_locus)
                tail_loci.append(tail_locus)
        return HomFam(head,head_loci,copy_n,new_comment),HomFam(tail,tail_loci,copy_n,new_comment)
    # Deletes loci of length below a given threshold
    def filter(self,threshold=-1):
        if threshold > 0:
            short_loci_list = []
            for i in range(len(self.loci)):
                if abs(self.loci[i].end-self.loci[i].start+1) < threshold:
                    short_loci_list.append(i)
            for i in sorted(short_loci_list,reverse=True):
                self.loci.pop(i)
    # compares two HomFam
    # other - HomFam
    # Return - int: compare value
    def __cmp__(self, other):
        if type(other) is HomFam:
            return cmp(self.id, other.id)
        else:
            return 1
        #endif
    #endfor
    
#    # tests equality of HomFam and string
#    # other - string
#    # Return - boolean: True if other = self.id and False otherwise
#    def __eq__(self, other):
#        if type(other) is str:
#            return self.id == other
#        else:
#            return self == other
#        #endif
#    #enddef
#
    # Redefinition of __hash__ so that class object remains hashable
#    def __hash__(self):
#        return hash((self.id,frozenset(self.loci)))
#    #enddef

    # returns the opening string of the HomFam in the form
    # \><id> <copy number> #<comment>
    # Return - str
    def opening_string(self):
        if len(self.comment) > 0:
            comment_str = " #" + self.comment
        else:
            comment_str = ""
        #endif
        
        if self.copy_number != 1:
            copy_str = " " +  str(self.copy_number)
        else:
            copy_str = ''
        #endif
        
        return ">" + self.id + copy_str + comment_str
    #enddef
    
    # returns a string representing itself in the format
    # <id> <copy number> #<comment>
    # locus1
    # locus2
    # ...
    # locusn
    # Return - str
    def __str__(self):
        string = self.opening_string() + "\n"
        
        for l in self.loci:
            string = string + str(l) + "\n"
        #endfor
        
        return string
    #enddef
    
    # writes the string representation to file
    # see __str__ for format
    # file_stream - file: file stream to write to
    def to_file(self, file_stream):
        #assert type(file_stream) is file
    
        file_stream.write(self.opening_string() + "\n")
            
        for l in self.loci:
            file_stream.write(str(l) + "\n")
        #endfor
    #enddef
    
    # reads the HomFam from a file
    # file_stream - the file_stream to read from
    # first_line - the first line of the file_stream if already read from otherwise None
    # Return - HomFam
    @staticmethod
    def from_file(file_stream, first_line,threshold=-1):
        #assert type(file_stream) is file
        #assert type(first_line) is str
    
        if first_line == None:
            line = file_stream.readline()
        else:
            line = first_line
        #endif
        
        if len(line) <= 0:
            return None, line
        #endif
        
        hom_fam = HomFam.opening_from_string(line)
                
        line = file_stream.readline()
        
        while len(line) > 0:
            line = line.strip()
            
            if len(line) > 0:
                if line[0]== '>':
                    break
                #endif
                
                if hom_fam != None:
                    locus = Locus.from_string(line,threshold)
                    
                    if locus != None:       
                        hom_fam.loci.append(locus)
                    #endif
                #endif
            #endif
                
            line = file_stream.readline()
        #endwhile
        
        return hom_fam, line
    #enddef
    
    # parses a HomFam from a string
    # see opening_str for format
    # string - str: the string to parse
    @staticmethod
    def opening_from_string(string):
        trunc_line = string.strip()
    
        if len(trunc_line) < 0 or trunc_line[0] != '>':
            print("Warning: not a hom. family. String: \'" + string + "\'")
            
            return None
        #endfor
            
        comment_split = trunc_line.split('#', 1)
        if len(comment_split) >= 2:
            comment = comment_split[1]
        else:
            comment = ""
        #endif
        trunc_line = comment_split[0].strip()
        
        split = trunc_line.split(' ', 1)
        
        if len(split[0]) > 0:
            ident = split[0][1:]
        else:
            print("Warning: hom. family missing identity: \'" + string + "\'")
            
            return None
        #endif
        
        copy_number = 1
        if len(split) > 1:
            if len(split[1]) > 0:
                copy_number = split[1].strip()
                                
                try:
                    copy_number = int(copy_number)
                except:
                    print("Warning: copy number not an integer for homology family. String : \'" + string + "\'")
                    return None
                #endtry
            else:
                print("Warning: unknown information in hom. family. String: \'" + string + "\'")
            #endif
        #endif
        
        return HomFam(ident, [], copy_number, comment)
    #enddef
#endclass

# reads hom. families from a file
# file_name - str: the name of the file to read from
# hom_fam_list - list of HomFam: the list to add to (Default = [])
# Return - list of HomFam: the list of hom. familes read
def read_hom_families_file(file_name,threshold=-1):
    file_stream = open(file_name)
    
    line = file_stream.readline()
    hom_fam_list = {}
#    hom_fam_list = []
    while len(line) > 0:
        trunc_line = line.strip()
        
        if len(trunc_line) > 0:
            if trunc_line[0] == '>':
                # read first hom. family
                hom_fam, line = HomFam.from_file(file_stream, trunc_line,threshold)
                if hom_fam != None:
                    hom_fam_list[hom_fam.id] = hom_fam
#                    hom_fam_list.append(hom_fam)
                #endif
                            
                # read the rest of the hom. families
                while len(line) > 0:
                    hom_fam, line = HomFam.from_file(file_stream, line,threshold)
                    if hom_fam != None and len(hom_fam.loci) > 0:
                        hom_fam_list[hom_fam.id] = hom_fam
#                        hom_fam_list.append(hom_fam)
                    #endif
                #endwhile
            else:
                line = file_stream.readline()
            #endif
        else:
            line = file_stream.readline()
        #endif
    #endif
    
    file_stream.close()
    
    return hom_fam_list
#enddef

# Double homologous families: take in a dictionary of 
# homologous families, and return heads and tails for 
# each locus of each family, identified with the 
# head/tail of said family. If the marker have no 
# orientation, exit with an error.

def double_hom_families(hom_families,num=0):
    doubled_families = {}
    identifiers = {}
    for b in hom_families:
        h,t = hom_families[b].doubler(num)
        identifiers[b] = (h.id,t.id)
        doubled_families[h.id] = h
        doubled_families[t.id] = t
    return doubled_families,identifiers

# writes hom. families to a file
# file_name - str: the name of the file to write to
# hom_fam_list - list of HomFam: the list to write (Default = [])
def write_hom_families_file(file_name, hom_fam_list):   
    file_stream = open(file_name, 'w')
    for hom_fam in hom_fam_list:
        try:
            file_stream.write(str(hom_fam_list[hom_fam]) + '\n')
        except:
            file_stream.write(str(hom_fam) + '\n')
        file_stream.flush()
    #endif
    file_stream.close()
#enddef

# finds a HomFam in a list given a marker id
# hom_fam_list - list of HomFam: the list to look in
# marker_id - str: the marker_id to look for
# Return - HomFam: the hom. family with marker id, marker_id
def find_marker(hom_fam_list, marker_id):
    try:
        index = hom_fam_list.index(marker_id)
    
        return hom_fam_list[index]
    except:
        return None
    #endtry 
#enddef
