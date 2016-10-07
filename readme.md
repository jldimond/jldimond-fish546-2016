#FISH546 2016 project repository

This is a repository for my 2016 FISH 546 project, involving analysis of _Porites_ spp. ddRAD-seq and EpiRAD-seq data.

I have recently obtained [ddRAD-seq](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037135) and [EpiRAD-seq](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12435/abstract) data for 48 samples of _Porites_ spp. corals collected in Belize. Most of these samples consist of branching _Porites_ spp., for which there is taxonomic uncertainty whether these comprise 3 different species or a single polymorphic species. My goal is to evaluate (1) if there is any genetic structuring of these individuals using the ddRAD data and (2) if there is any epigenetic structuring of these individuals using the EpiRAD data. I will be using iPyrad to obtain the clusters of loci, and will be filtering out symbiont sequences by using the "denovo-reference" option to obtain sequences that do not match the Symbiodinium minutum genome.

The directory structure is as follows:

'analyses'  - Files resulting from analyses.
'data' -  Information and links to methods and raw data.
'notebooks'	- Jupyter notebooks.
'scripts' - Scripts such as R scripts used for analyses.