# CIMICE pipeline 

This is the CIMICE pipeline for bacterial genome analysis that is associated with th tool 
CIMICE ('https://github.com/redsnic/tumorEvolutionWithMarkovChains').
It's goal is to provide an easy to use software to extract mutational
matrices from annotated bacterial sequences. 

Snippy ('https://github.com/tseemann/snippy') is used for variant calling 
and is needed for the full analysis, the version I used is that from anaconda2. 

## Installation

This software is completely written in python and needs no installazion, you 
can download the sources directly from ('https://github.com/redsnic/CIMICEPipeline/') or,
if you have git installed, issuing the command 

```
git clone https://github.com/redsnic/CIMICEPipeline/
```

## Usage

To run the pipeline write and execute the following in a terminal 

```
python3 <path_to_CIMICEpipeline.py> <args>
```

In this table you can see all the supported arguments and how they are formatted

```
========================================================================================================

USAGE: CIMICEpipeline <MODE> <ARGS>

========================================================================================================

MODES: 

   FULL [CLEAN] <path_to_reference (.gbk)> <path_to_sequences (.fa)> 
                [WITHOUT] [mutation_kind1 mutation_kind2 ...]        
   POSTPROCESS [ONLY_GENES] <reference_identifier> <path_to_csv_folder>
                [WITHOUT] [mutation_kind1 mutation_kind2 ...]        
   SORT <dataset to be sorted>
   DIFF id <references_folder> <samples_folder>
   MERGE <dataset1 dataset2 ...>
   HELP

========================================================================================================
```

### Supported execution options: 

   #### FULL: 
         execute snippy for variant calling on a reference and a set of sequences. It is possible
         to choose the kinds of mutation to be used in the post processing (or not to if
         the without keyword is used. The output are a dataset in TRONCO format and a folder
         containing the csv files computed by snippy (placed in the current directory).
         CLEAN option can be used to erase unnecessary snippy output
         Note: the execution of this mode could take a long time; for reference, to analyze 
         E.coli genomes (~ 5.4 Mbp) on an four core machine, the time necessary 
         for the execution could be estimated with this function:
         time = (2 +/- 1 minutes) x (number_of_references) x (number_of_samples)

   #### POSTPROCESS: 
         run only the post processing activities. It takes in input an identifier and 
         the folder containing the output csv files from snippy. It is possible 
         to choose the kinds of mutation to be used in the post processing (or not to if
         the without keyword is used. The output is a dataset in TRONCO format.
         If ONLY_GENES option is enabled LOCUS_TAG information will not be used if 
         gene name is not present.

   #### SORT: 
         sort a TRONCO formatted dataset so that genes and samples are sorted by most mutates.

   #### DIFF: 
         prepare difference datasets using gene information from .gbk/.gbff files.
         Output is presented in two files: dataset_id_absent and dataset_id_present,
         they are in TRONCO format and contain informations about the gene gains and 
         and loss in the given samples if compared to the given references.  

   #### MERGE: 
        merge the informations of a list of datasets in a single one called 'merged_datasets.CAPRI'

   #### HELP:
        print this message on the standard output.

### Notes

The arguments in '<>' parenthesis are mandatory while that in '[]' are optional.
If no mutation_kind is given the default setting 'WITHOUT synonymous_variant' is used.
To take into account every kind of mutations use 'WITHOUT' keyword without specifying any
mutation_kind. Refer to VCF format ANN tag for a full list of mutation types, most common are: 
synonymous_variant, missense_variant, frameshift_variant and stop_gained.
To analyze the output dataset in CIMICE remember to use -c or --capri options.
