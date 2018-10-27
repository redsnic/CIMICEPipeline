'''
Created on 25 ott 2018

@author: redsnic
'''

import fileinput,os


def transpose(matrix):
    return [[matrix[i][j] for i in range(len(matrix))] for j in range(len(matrix[0]))]

class Dataset(object):
    '''
    This class is meant to model CIMICE/TRONCO input files and to 
    do simple manipulation on them
    '''

    def toCAPRI(self):
        '''
        returns a string of this dataset shown in CAPRI format
        '''
        ret = ""
        ret += "s\\g" + " " + " ".join(self.genes) + "\n"
        for i in range(len(self.samples)):
            ret += os.path.basename(self.samples[i]) + " " + " ".join(map(str,self.matrix[i])) + "\n"
        return ret
        
    
    def _sort_by_sum(self, labels):
        '''
        sorting lines of a matrix by their sum
        labels: identifier of the lines to be kept in the correct position  
        '''
        tup = []
        for i in range(len(labels)):
            tup.append((labels[i], sum(self.matrix[i]), self.matrix[i]))
        tup.sort(key = lambda t:t[1], reverse=True)
        
        for i in range(len(labels)):
            labels[i] = tup[i][0]
            self.matrix[i] = tup[i][2]
            
    def getDictionaries(self):
        out = []
        for i in range(len(self.samples)):
            dic = {}
            for j in range(len(self.genes)):
                if self.matrix[i][j] == 1:
                    dic[self.genes[j]] = self.matrix[i][j]
            out.append((self.samples[i], dic))
        return out
        
        
    def sort(self):
        '''
        sort this dataset to:
        order samples from most mutated to less mutated and
        order genes from most mutated to less mutated
        (useful for visualization)
        '''
        if len(self.genes) > 0:
            self._sort_by_sum(self.samples)
            self.matrix = transpose(self.matrix)
            self._sort_by_sum(self.genes)
            self.matrix = transpose(self.matrix)
        else:
            print("Warning: empty dataset")
    def read(self, path):
        '''
        read input dataset from path
        '''
        first=True
        for line in fileinput.input(path):
            if(first):
                self.genes = line.split()[1:]
                first=False
            else:
                l = line.split()
                self.samples.append(l[0])
                self.matrix.append(list(map(int, l[1:])))
                
    def create(self, genes, lines):
        '''
        initialize dataset 
        genes : list of genes
        lines : sample, genes as 0/1 (absence/presence) 
        '''
        self.genes = genes
        for x in lines:
            self.samples.append(x[0])
            self.matrix.append(x[1:])
        

    def __init__(self):
        '''
        Constructor
        '''
        self.genes = []
        self.samples = []
        self.matrix = []
        
    
if __name__ == "__main__":
    data = Dataset()
    data.read("/home/redsnic/CIMICE/tumorEvolutionWithMarkovChains/datasets/dataset_NC_018658.CAPRI")    
    print(data.toCAPRI())
    data.sort()
    print(data.toCAPRI())
    
    
    