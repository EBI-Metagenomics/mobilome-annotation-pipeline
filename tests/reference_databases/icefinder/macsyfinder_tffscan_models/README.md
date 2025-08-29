# TFFscan: MacSyFinder models for the annotation of the type IV filament super-family

This set of MacSyFinder's models are dedicated to the genomic detection of the cellular machineries that are members of the type IV filament super-family. 
The systems that can be detected are the following: 

- Archaeal-T4P: the set of archeal type IV pili, comprising the archaeal flagellum (or archaeallum), the UV-induced pilus Ups, the bindosome Bas, etc... 
- ComM: the competence machinery found in monoderms, including the classical Bacillus machinery
- MSH: the mannose-sensitive hemagglutinin pilus
- T2SS: the type II secretion system 
- T4aP: the type IV pilus, sometimes also called the Type IVa pilus
- T4bP: the type IVb pilus, comprising the R64 thin pilus, toxin-coregulated pilus (TCP), bundle-forming pilus (Bfp), longus pilus, and Cof pilus
- Tad: the Tad (thigh-adherence) pilus


This set of models corresponds to those published in 2019 in Plos Biology (see below for full reference), updated to the format to enable systems search with **MacSyFinder version 2**.

## Installation and Usage with MacSyFinder

First, the `macsyfinder` program should be [installed](http://macsyfinder.readthedocs.io/en/latest/). This will also install the `macsydata` tool that enables to easily install this package of MacSyFinder models from this repository. 


The basic commands to run are then:

    macsydata install TFFscan


to install the TFFscan package. 

    macsyfinder --db-type ordered_replicon \
		--sequence-db myproteins.fasta \
		--models TFFscan system 		


to run the search on your favorite organism's genom, where `system` is one or multiple systems listed above, or `all` to search for all the above listed systems
(see [MacSyFinder's documentation](http://macsyfinder.readthedocs.io/en/latest/)). 


It has to be noted that to ensure the highest annotation specificity, it is recommended to search for **all** systems in the package at once. 


## References

- Abby Sophie S., Néron Bertrand, Ménager Hervé, Touchon Marie, Rocha Eduardo P. C.
  (2014).
  MacSyFinder: A Program to Mine Genomes for Molecular Systems with an Application to CRISPR-Cas Systems.
  PLoS ONE, 9 (10), pp. e110726.
  [doi:10.1371/journal.pone.0110726](http://dx.doi.org/10.1371/journal.pone.0110726)


- Denise Rémi, Abby Sophie S., Rocha Eduardo P. C. (2019). 
  Diversification of the type IV filament superfamily into machines for adhesion, protein secretion, DNA uptake, and motility.
  PLoS Biology 17(7): e3000390.
  [doi:10.1371/journal.pbio.3000390](https://doi.org/10.1371/journal.pbio.3000390)
  

- Abby Sophie S., Cury Jean, Guglielmini Julien, Néron Bertrand, Touchon Marie, Rocha Eduardo P. C.
  (2016).
  Identification of protein secretion systems in bacterial genomes.
  Scientific Reports, 6, pp. 23080.
  [doi:10.1038/srep23080](http://dx.doi.org/10.1038/srep23080)

For more details on this set of models, see the following reference: [Denise et al. 2019](https://doi.org/10.1371/journal.pbio.3000390), Diversification of the type IV filament superfamily into machines for adhesion, protein secretion, DNA uptake, and motility, PLoS Biology. 

NB: The HMM protein profiles included were either built de novo or obtained from a previous publication. Check the above references and metadata file for more details.
