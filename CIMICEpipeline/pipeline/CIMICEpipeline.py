'''
Created on 25 ott 2018

@author: redsnic

This file contains the I/O and the settings management procedures of the CIMICE pipeline

'''

import Dataset, os, shutil, subprocess, findDifferences
import SnippyCsv
from mergeDatasets import mergeDatasets

def help_msg():
    msg = "\n"
    msg += "===================================         CIMICE pipeline           ===================================\n"
    msg += "\n"
    msg += "This is the CIMICE pipeline for bacterial genome analysis.\n"
    msg += "It's goal is to provide an easy to use software to extract mutational\n"
    msg += "matrices from annotated bacterial sequences. \n"
    msg += "Snippy ('https://github.com/tseemann/snippy') is used for variant calling \n"
    msg += "and is needed for the full analysis, the version I used is that from anaconda2. \n"
    msg += "\n"
    msg += "========================================================================================================\n"
    msg += "\n"
    msg += "USAGE: CIMICEpipeline <MODE> <ARGS>\n"
    msg += "\n"
    msg += "========================================================================================================\n"
    msg += "\n"
    msg += "MODES: \n\n"
    msg += "   FULL [CLEAN] <path_to_reference (.gbk)> <path_to_sequences (.fa)> \n"
    msg += "                [WITHOUT] [mutation_kind1 mutation_kind2 ...]        \n"
    msg += "   POSTPROCESS [ONLY_GENES] <reference_identifier> <path_to_csv_folder>\n"
    msg += "                [WITHOUT] [mutation_kind1 mutation_kind2 ...]        \n"
    msg += "   SORT <dataset to be sorted>\n"
    msg += "   DIFF id <references_folder> <samples_folder>\n"
    msg += "   MERGE <dataset1 dataset2 ...>\n"
    msg += "   HELP\n"
    msg += "\n"
    msg += "========================================================================================================\n"
    msg += "\n"
    msg += "DESCRIPTIONS: \n\n"
    msg += "   FULL: execute snippy for variant calling on a reference and a set of sequences. It is possible\n"
    msg += "         to choose the kinds of mutation to be used in the post processing (or not to if\n"
    msg += "         the without keyword is used. The output are a dataset in TRONCO format and a folder\n"
    msg += "         containing the csv files computed by snippy (placed in the current directory).\n"
    msg += "         CLEAN option can be used to erase unnecessary snippy output\n"
    msg += "         Note: the execution of this mode could take a long time; for reference, to analyze \n"
    msg += "         E.coli genomes (~ 5.4 Mbp) on an four core machine, the time necessary \n"
    msg += "         for the execution could be estimated with this function:\n"
    msg += "         time = (2 +/- 1 minutes) x (number_of_references) x (number_of_samples)\n"
    msg += "\n"
    msg += "   POSTPROCESS: run only the post processing activities. It takes in input an identifier and \n"
    msg += "         the folder containing the output csv files from snippy. It is possible \n"
    msg += "         to choose the kinds of mutation to be used in the post processing (or not to if\n"
    msg += "         the without keyword is used. The output is a dataset in TRONCO format.\n"
    msg += "         If ONLY_GENES option is enabled LOCUS_TAG information will not be used if \n"
    msg += "         gene name is not present.\n"
    msg += "\n"
    msg += "   SORT: sort a TRONCO formatted dataset so that genes and samples are sorted by most mutates.\n"
    msg += "\n"
    msg += "   DIFF: prepare difference datasets using gene information from .gbk/.gbff files.\n"
    msg += "         Output is presented in two files: dataset_id_absent and dataset_id_present,\n"
    msg += "         they are in TRONCO format and contain informations about the gene gains and \n"
    msg += "         and loss in the given samples if compared to the given references.  \n"
    msg += "\n"
    msg += "   MERGE: merge the informations of a list of datasets in a single one called 'merged_datasets.CAPRI'\n"
    msg += "\n"
    msg += "   HELP: print this message.\n"
    msg += "\n"
    msg += "Note that arguments in '<>' parenthesis are mandatory while that in '[]' are optional.\n"
    msg += "If no mutation_kind is given the default setting 'WITHOUT synonymous_variant' is used.\n"
    msg += "To take into account every kind of mutations use 'WITHOUT' keyword without specifying any\n"
    msg += "mutation_kind. Refer to VCF format ANN tag for a full list of mutation types, most common are: \n"
    msg += "synonymous_variant, missense_variant, frameshift_variant and stop_gained.\n"
    msg += "To analyze the output dataset in CIMICE remember to use -c or --capri options.\n"
    msg += "\n"
    msg += "========================================================================================================\n"
    
    return msg

def to_dataset(out):
    '''
    converts a list (sample, genes ...) into a dataset
    '''
    all_genes = []
    samples = []
    dicts = []
    #prepare lines
    for line in out:
        samples.append(line[0])
        all_genes += line[1:]
        dic = {}
        for x in line:
            dic[x] = 1
        dicts.append(dic)
    # select unique genes
    all_genes = list(set(all_genes))
    # prepare lines
    lines = []
    for i in range(len(samples)):
        line = [samples[i]]
        for gene in all_genes:
            try:
                line.append(dicts[i][gene])
            except KeyError:
                line.append(0)
        lines.append(line)
        
    data = Dataset.Dataset()   
    data.create(all_genes,lines) 
        
    return data 
        

class CIMICEpipline():
    '''
    Class to manage the CIMICE pipeline execution
    '''
    
    def full(self, args):
        '''
        run snippy and the post processing 
        reference_path data_path 
        '''
        
        #cleanup flag TODO
        self.remove_temp_files = False
        if len(args) > 0 and args[0].upper() == "CLEAN":
            self.remove_temp_files = True
            args = args[1:]
            
        
        if len(args) < 2:
            print("Error: wrong number of arguments for snippy execution")
            sys.exit(-3)
            
        
        print("Preparing snippy execution...")
        
        folder = ""
        reference_name = os.path.splitext(os.path.basename(args[0]))[0]
        
        try:
            
            #create output folder for csvs
            folder = 'csvs_' + reference_name
            
            if os.path.exists(folder):
                shutil.rmtree(folder)
            os.makedirs(folder) 
            
        except:
            print("Error in the creation of csv files folder")
            sys.exit(-3)
                
        
        output_dir_prefix = reference_name+"_sinppy"
        try: 
            for seq in os.listdir(args[1]):
                print("Analyzing " + seq + " ...")
                extension = os.path.splitext(os.path.basename(seq))[-1]
                if extension == ".fa" or extension == ".fasta" or extension == ".fna":
                    subprocess.call(["snippy","--outdir", output_dir_prefix+"_"+seq+"_elaboration", "--ref", args[0], "--ctgs", args[1]+ "/" + seq , "--quiet", "--force"])
                    shutil.move(output_dir_prefix+"_"+seq+"_elaboration" + "/" + "snps.csv", folder + "/" + os.path.splitext(os.path.basename(seq))[0] + ".csv" )
                    if os.path.exists(output_dir_prefix+"_"+seq+"_elaboration") and self.remove_temp_files:
                        shutil.rmtree(output_dir_prefix+"_"+seq+"_elaboration")
        except:
            print("Error in snippy execution")
            sys.exit(-4)
        
         
        if len(args[2:])>0:
            self.post_process([reference_name, folder] + args[2:])
        else:
            self.post_process([reference_name, folder, "WITHOUT", "synonymous_variant"])
            
      
        print("done, procedure completed successfully")
    
    def sort(self, args):
        '''
        sorting procedure
        '''
        print("sorting phase:")
        if len(args) != 1:
            print("Error: wrong arguments for sorting mode")
            sys.exit(-1)
        try:        
            data = Dataset.Dataset()
            data.read(args[0])
            data.sort()
            f = open(args[0], "w")
            f.write(data.toCAPRI())
        except:
            print("Error reading the DATASET, check your path")
            sys.exit(-2)
        
    def post_process(self, args):
        
        if(len(args) < 2):
            print("Error: too few argument for post processing")
            sys.exit(-2)
        
        use_locus_tag_if_necessary = True
        if args[0].upper() == "ONLY_GENES":
            use_locus_tag_if_necessary = False
            args = args[1:]
         
        ref_id = args[0]    
            
        print("post processing phase:")
        # prepare output folder
        folder = 'tables_'+ref_id
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.makedirs(folder)  
        out = []
        try:        
            file_list = os.listdir(args[1])
            for el in file_list:
                if el.endswith(".csv"):
                    print("Elaborating CSV: " + str(el))
                    csv = SnippyCsv.SnippyCsv(args[1] + "/" + el, args[2:], use_locus_tag_if_necessary)
                    out.append(csv.to_line())
        except:
            print("Error reading the csv files, check your path")
            sys.exit(-2)
        
        # prepare for output
        f = open("dataset_"+args[0]+".CAPRI", "w")
        dataset = to_dataset(out)
        dataset.sort()
        f.write(dataset.toCAPRI())            
        
        
    
    def __init__(self):
        '''
        constructor
        '''
        pass

import sys

if __name__ == "__main__":
    
    execution_type = ""
    
    try:
        execution_type = sys.argv[1]
    except IndexError:
        print(help_msg()) 
        sys.exit(-1)
        
    execution_type = execution_type.upper()
    
    pipeline = CIMICEpipline()    
    
    if execution_type == "FULL" or execution_type == "F":
        pipeline.full(sys.argv[2:])
    elif execution_type == "POST" or execution_type == "P" or execution_type == "POSTPROCESS":
        pipeline.post_process(sys.argv[2:])
    elif execution_type == "SNIPPY":
        pipeline.snippy(sys.argv[2:])
    elif execution_type == "SORT":
        pipeline.sort(sys.argv[2:])
    elif execution_type == "DIFF":
        if(len(sys.argv) == 5):
            findDifferences(sys.argv[3], sys.argv[4], sys.argv[2])
        else:
            print("Not enough arguments for DIFF mode")
            sys.exit(-1)
    elif execution_type == "MERGE":
        if(len(sys.argv) > 3):
            mergeDatasets(sys.argv[2:])
        else:
            print("Not enough arguments for MERGE mode")
            sys.exit(-1)
    elif execution_type == "HELP": 
        print(help_msg())  
    else:
        print("No valid option specified, run 'python3 CIMICEpipeline help' for the list of valid options and arguments")
    
    
    