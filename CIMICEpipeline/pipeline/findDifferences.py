'''
Created on 27 ott 2018

@author: redsnic

difference dataset preparation
'''

import os, Dataset, fileinput

def extract_genes(file_gbk):
    l = []
    print("reading file "+file_gbk)
    for line in fileinput.input(file_gbk):        
        if "/gene=" in line:
            l.append(line.split()[0][7:].replace('"', ''))
    return set(l)        

def findDifferences(ref_folder, seq_folder, identifier):
    '''
    creates two dataset containing the set of genes absent in 
    or present in a set of references and sequences
    '''
    print("preparing difference datasets...")
    all_genes = set()
    ref_genes = set()
    for ref in os.listdir(ref_folder):
        ref_genes.update(extract_genes(ref_folder + "/" + ref))
    
    all_genes.update(ref_genes)
    
    samples = []
    genes_tables = []
    genes_sets = []
    for seq in os.listdir(seq_folder):
        samples.append(seq)
        genes = extract_genes(seq_folder + "/" +seq)
        all_genes.update(genes)
        genes_sets.append(genes)
        genes_table = {}
        for g in genes:
            genes_table[g] = 1 
        genes_tables.append(genes_table)
    
    
    absent = Dataset.Dataset()
    present = Dataset.Dataset()
    
    i = 0
    present_lines = []
    absent_lines = []
    for genes_set in genes_sets:
        absent_list = ref_genes.difference(genes_set)
        
        absent_table = {}
        for a in absent_list:
            absent_table[a] = 1 
        
        absent_line = [os.path.splitext(samples[i])[0]]
        for g in all_genes:
            try:
                absent_line.append(absent_table[g])
            except KeyError:
                absent_line.append(0)
        absent_lines.append(absent_line)
        
        present_list = genes_set.difference(ref_genes)
        
        present_table = {}
        for a in present_list:
            present_table[a] = 1
        
        present_line = [os.path.splitext(samples[i])[0]]
        for g in all_genes:
            try:
                present_line.append(present_table[g])
            except KeyError:
                present_line.append(0)
        present_lines.append(present_line)
        
        
        i+=1
    
         
    absent.create(list(all_genes), absent_lines)
    absent.sort()
    f = open("dataset_"+identifier+"_absent.CAPRI", "w")
    f.write(absent.toCAPRI())
    
    present.create(list(all_genes), present_lines)  
    present.sort()  
    f = open("dataset_"+identifier+"_present.CAPRI", "w")
    f.write(present.toCAPRI())    
    print("done, output in files " + "dataset_"+identifier+"_absent.CAPRI" + " and " + "dataset_"+identifier+"_present.CAPRI")

if __name__ == "__main__":
    refs = "/home/redsnic/CIMICEPipeline/PRJNA285020/references"
    seqs = "/home/redsnic/CIMICEPipeline/PRJNA285020/GBK"
    findDifferences(refs, seqs, "TEST")


