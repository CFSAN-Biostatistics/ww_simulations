# ww_simulations

A series of scripts to generate simulated datasets of different variant composition for SC2 wastewater surveillance efforts.

This is in active development and not ready for primetime yet...hopefully soon.

# Dependencies
* python (3.8.1)
* [art](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) 
* [in_silico_pcr.pl](https://github.com/egonozer/in_silico_pcr)
    
# Usage
* `generate_replicates.py` can be used to create replicate datasets each with a random composition of a fixed number of variants

Example: `generate_replicates.py -i epi_isls.txt -n 2 -o test_output`

* `generate_simulated_datasets.py` takes the ouput of `generate_replicates.py` and a number of other files/arguments and generates simulated data reflecting the composition of the variants
 
Example: ```generate_simulated_datasets.py -f relevantGISAID.fa -i test_output_replicate_1.tsv -p neb_vss1a_primer.tsv -n 10```

# Workflow/strategy
* `generate_simulated_datasets.py` models differences in variation composition by multiplying the number of each amplicon from each variant by the percent abundance it is in the simulated 'sample'. 
* The total reads of is determined by the `-f` flag passed to `art` which is "the fold of read coverage to be simulated or number of reads/read pairs generated for each amplicon"


# License
See the LICENSE.txt file included in the SNP Pipeline distribution.
