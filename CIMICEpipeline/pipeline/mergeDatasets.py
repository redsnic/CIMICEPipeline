'''
Created on 27 ott 2018

@author: redsnic

dataset merging procedure
'''

import Dataset


def mergeDatasets(paths):
    print("Merging datasets...")
    datasets = []
    for path in paths:
        d = Dataset.Dataset()
        d.read(path)
        datasets.append(d)
    
    total_genes = set()
    samples_info = {}
    for d in datasets:
        total_genes.update(d.genes)
        for s in d.samples:
            samples_info[s] = {}
    
    total_samples = set()
    for d in datasets:
        dicts = d.getDictionaries()
        for p in dicts:
            sample = p[0]
            total_samples.add(sample)
            #second element is always one
            for x,_ in p[1].items():
                try:
                    samples_info[sample][x] = 1
                except KeyError:
                    samples_info[sample] = {}
                    samples_info[sample][x] = 1
    

    lines = []
    for s in total_samples:
        line = [s]
        for g in total_genes:
            try:
                line.append(samples_info[s][g])
            except KeyError:
                line.append(0)
        lines.append(line)
        
    final_dataset = Dataset.Dataset()
    
    final_dataset.create(list(total_genes), lines)   
    final_dataset.sort()
    
    f = open("merged_datasets.CAPRI", "w")
    f.write(final_dataset.toCAPRI())  
    print("done, output in file merged_datasets.CAPRI")  
    
if __name__ == "__main__":
    mergeDatasets(["/home/redsnic/eclipse-workspace/CIMICEpipeline/dataset_test_absent.CAPRI", "/home/redsnic/eclipse-workspace/CIMICEpipeline/dataset_test_present.CAPRI", "/home/redsnic/eclipse-workspace/CIMICEpipeline/dataset_CP008957.CAPRI"])
         
    