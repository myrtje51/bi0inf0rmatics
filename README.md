# bi0inf0rmatics

How many functions does the program contain? 
Functions for mapping DrugBank and STITCH: 
- STITCHfilter() --> Filters STITCH dataset bij a confidence score of 0.7 or higher. 
- STITCHgetInchi() --> Get the InchiKey's using PubChem (input: CID, output: InchiKey) (Hopefully we'll be able to use the package pubchempy for that).
- STITCHgetDBID() --> Get the DrugBank ID's using Unichem (input: InchiKey, output: DrugBank ID) (there isn't a package on this). 
- MapDBandSTITCH() --> Map the DrugBank dataset to the STITCH dataset using the DrugBank ID (R could be used for this). 

Funtions for the Ageing Clusters Resource:
- ACRfilter() --> Filters the ACR by only keeping the genes that are present in two or more of the categories (you can use some if-statements for that probably since it's not a big dataset). 
- ACR_DTenrich() --> The filtered dataset is enriched using the Drug & Target dataset. This can be done using R with the Fisher's exact test. 

Functions for the different biological levels:
- Reactome_ACRenrich() --> The filtered geneset and the reactome terms are enriched against eachother. For this a R-package is used called: EnrichPathway.
- Reactome_DTenrich() --> 

- KEGG_ACRenrich() --> The filtered geneset and the KEGG terms are enriched against eachother. For this a R-package is used as well called: EnrichKEGG. 
- KEGG_DTenrich() --> 

- GO_ACRenrich() --> The filtered geneset and the GO terms are enriched against eachother. For this a R-package is used called: EnrichGO. 
- GO_DTenrich() --> 

- PPI_ACRenrich() --> 
- PPI_DTenrich() --> 

- FirstRanking() --> Out of each enrichment comes a list of drugs with p-values. These p-values will determine the ranking of the drugs. The smaller the p-value the better (so a drug with a small p-value will be high up in the ranking).  

Functions for getting the last ranking: 
- Last_ranking() --> 
