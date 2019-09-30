# bi0inf0rmatics

How many functions does the program contain? 

Functions for mapping DrugBank and STITCH: 
- STITCHfilter() --> Filters STITCH dataset by a confidence score of 0.7 or higher. This could be done using import sys and then taking the last column. If the value in that column is higher than or equal to 700 the line needs to be written to another file.  
- STITCHgetInchi() --> Get the InchiKey's using PubChem (input: CID, output: InchiKey) (Hopefully we'll be able to use the package pubchempy for that).
- STITCHgetDBID() --> Get the DrugBank ID's using Unichem (input: InchiKey, output: DrugBank ID) (there isn't a package on this). Maybe you don't have to use Unichem because DrugBank just made a new dataset that has the DrugBank ID and the InchiKey so maybe we could map those two datasets together? So if the InchiKey is a certain "value" then write DrugBank ID to another file (or add to certain table).  
- MapDBandSTITCH() --> Map the DrugBank dataset to the STITCH dataset using the DrugBank ID (R could be used for this). So you can use %in% for that probably.  

Funtions for the Ageing Clusters Resource:
- ACRfilter() --> Filters the ACR by only keeping the genes that are present in two or more of the categories (you can use some if-statements for that probably since it's not a big dataset). There is a column that counts in how many categories the gene is present so you can just tell the program that if the value in that column is higher than 1 write that to another file.  
- ACR_DTenrich() --> The filtered dataset is enriched using the Drug & Target dataset. This can be done using R with the Fisher's exact test. 

Functions for PPI filter: 
- PPIfilter() --> Filters PPI dataset by a confidence score of 0.9 or higher. This could be done using import sys en then taking the last column. If the value in that column is higher than or equal to 900 the line needs to be written to another file. 

Functions for the different biological levels:
- Reactome_ACRenrich() --> The filtered geneset and the reactome terms are enriched against eachother. For this a R-package is used called: EnrichPathway. The output will be a list of age-related reactome terms. 
- Reactome_DTenrich() --> The list of age-related reactome terms gets enriched against the Drug & Target dataset. There is no R-package for that so that just needs to be done with the fisher's exact test. The output will be a list of drugs. 

- KEGG_ACRenrich() --> The filtered geneset and the KEGG terms are enriched against eachother. For this a R-package is used as well called: EnrichKEGG. The output will be a list of age-related KEGG terms. 
- KEGG_DTenrich() --> The list of age-related KEGG terms gets enriched against the Drug & Target dataset. There is no R-package for that so that just needs to be done with the fisher's exact test. The output will be a list of drugs. 

- GO_ACRenrich() --> The filtered geneset and the GO terms are enriched against eachother. For this a R-package is used called: EnrichGO. The output will be a list of age-related GO terms.  
- GO_DTenrich() --> The list of age-related GO terms gets enriched against the Drug & Target dataset. There is no R-package for that so that just needs to be done with the fisher's exact test. The output will be a list of drugs. 

- PPI_ACRenrich() --> The filtered geneset and the PPI's are enriched against eachother. There is no R-package for this so only a fisher's exact test needs to be done. The output will be a list of age-related PPI's. 
- PPI_DTenrich() --> The list of age-related PPI's gets enriched against the Drug & Target dataset. There is no R-package for that so that just needs to be done with the fisher's exact test. The output will be a list of drugs. 

- FirstRanking() --> Out of each enrichment comes a list of drugs with p-values. These p-values will determine the ranking of the drugs. The smaller the p-value the better (so a drug with a small p-value will be high up in the ranking). Each list will be written to a file. So you will end up with 7 different files. This way, you will be able to take the ranking of each drug and take an average which will be important in the last function (CalLastRanking()).  

Functions for getting the last ranking: 
- CalLastRanking() --> From each drug a average ranking needs to be calculated. This ranking is the average of the ranking in each list (of the seven lists). You will first make a variable called ranking that starts at 0 and counts all the different rankings from one certain drug. Then you will tell the program to search for a certain drug by name in each file. After that you will tell the program to take the ranking and add the ranking to the rest of the rankings. Then you calculate the average and put that in a variable called average_ranking. This variable including the drug can then be written to a different file. 
- LastRanking() --> The list will be ranked by the average ranking. This could be done by the function series.rank() in pandas. 
