'''
Created on 25 ott 2018

@author: redsnic

Reader for snippy's csv files
'''

import fileinput, os

class SnippyCsv(object):
    '''
    Class to represent sippys' csv files with variant calling informations
    '''

    def __str__(self):
        ret = ""
        ret += ",".join(self.header) + "\n"
        for i in range(len(self.table["CHROM"])):
            ret += ",".join([self.table[h][i] for h in self.header])
            ret += "\n"
        return ret
    
    def to_line(self):
        '''
        returns a list: sample mutated genes 
        '''
        line = [self.name]
        alterations = []
        # fix bad choices for csv format
        for i in range(len(self.table["EFFECT"])):
            if self.table["EFFECT"][i] != "":
                self.table["EFFECT"][i] = self.table["EFFECT"][i].split()[0]
    
        # selection loop
        for i in range(len(self.table["CHROM"])):
            if self.accept_all or ((self.to_remove == False) and self.table["EFFECT"][i] in self.selected_kinds_of_mutation) or ((self.to_remove == True) and not(self.table["EFFECT"][i] in self.selected_kinds_of_mutation)):
                if self.table["GENE"][i] != "":
                    alterations.append(self.table["GENE"][i])
                elif self.use_locus_tag_if_necessary and self.table["LOCUS_TAG"][i] != "":
                    alterations.append(self.table["LOCUS_TAG"][i])
                # else skip mutation...
        # remove repetitions 
        line += list(set(alterations))
        return line

    def read(self, path):
        '''
        read csv from file (general for any csv)
        '''
        self.name = os.path.splitext(os.path.basename(path))[0]
        first=True
        self.table = {}
        self.header = []
        for line in fileinput.input(path):
            if(first):
                self.header = line.split(",")
                for el in self.header:
                    self.table[el] = []
                first=False
            else:
                l = line.split(",") 
                i=0
                for h in self.header:
                    try:
                        self.table[h].append(l[i].replace("\n", "")) # remove newline for last element
                        i+=1
                    except IndexError:
                        self.table[h].append("")
                        i+=1
        

    def __init__(self, path, selected_kinds_of_mutation = None, use_locus_tag_if_necessary = True):
        '''
        Constructor
        '''
        self.use_locus_tag_if_necessary = use_locus_tag_if_necessary
        self.accept_all = False
        self.to_remove = False
        if selected_kinds_of_mutation == None or len(selected_kinds_of_mutation) == 0:
            self.accept_all = True  # default
        elif selected_kinds_of_mutation[0].upper() == "WITHOUT":
            self.to_remove = True   # remove selected kinds of mutation
        
        self.selected_kinds_of_mutation = selected_kinds_of_mutation
        self.read(path)
        
        
    
        
        
        
        