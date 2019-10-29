# bi0inf0rmatics

For new users, some packages that are needed to run this program: 

Python packages: 
- rpy2
- pandas
- numpy
- scipy
- requests

R packages (within anaconda/miniconda): 
- clusterProfiler
- ReactomePA 
- STRINGdb
- org.Hs.eg.db

How many functions does the program contain? 

Functions for mapping DrugBank and STITCH: 
- STITCHfilter() --> Filters STITCH dataset by a confidence score of 0.7 or higher. This could be done using import sys and then taking the last column. If the value in that column is higher than or equal to 700 the line needs to be written to another file.  
- get_STITCH_inchikey --> Get the InchiKey's using PubChem (input: CID, output: InchiKey). We used the package requests to get and API through Python.  
- STITCHgetDBID() --> This function maps the STITCH database against the DrugBank database using the inchikey.

Funtions for the Ageing Clusters Resource:
- ACRfilter() --> Filters the ACR by only keeping the genes that are present in two or more of the categories (you can use some if-statements for that probably since it's not a big dataset). There is a column that counts in how many categories the gene is present so you can just tell the program that if the value in that column is higher than 1 write that to another file. 

Functions for PPI filter: 
- PPIfilter() --> Filters PPI dataset by a confidence score of 0.9 or higher. This could be done using import sys en then taking the last column. If the value in that column is higher than or equal to 900 the line needs to be written to another file. 

Functions for the conversion of several kinds of id's (in case you need them): 
- protein_to_entrez() --> converts the protein ensembl id's to entrez gene id's. 
- entrez_to_protein() --> converts the entrez gene id's to protein ensembl id's.

Functions for the full enrichment: 
- make_dictio_DT() --> makes a dictionary out of the dataset where DrugBank and STITCH are mapped. 
- cluster_profiler_KEGG() --> enriches the given gene list using the KEGG dataset. This is done with a R-package called: 
clusterProfiler. To be able to use this package in Python, rpy2 is used in this function. 
- cluster_profiler_GO_MF() --> enriches the given gene list using the GO molecular functions dataset. This is done with a R-
package called: clusterProfiler. To be able to use this package in Python, rpy2 is used in this function. 
- cluster_profiler_GO_CC() --> enriches the given gene list using the GO cellular component dataset. This is done with a R-
package called: clusterProfiler. To be able to use this package in Python, rpy2 is used in this function. 
- cluster_profiler_GO_BP() --> enriches the given gene list using the GO biological process dataset. This is done with a R-package called: clusterProfiler. To be able to use this package in Python, rpy2 is used in this function. 
- Reactome() --> enriches the given gene list using the Reactome dataset. This is done with a R-package called: ReactomePA.
To be able to use this package in Python, rpy2 is used in this function. 
- ppi_interactions() --> gets the protein-protein interactions from the R-package: STRINGdb. This package takes a gene list
and maps it against STRINGv10. 
- ppi_dictio() --> makes a dictionary out of the results that get out of the function: ppi_interactions(). 

Classes for the full enrichment:
- main_enrichments(object) --> does the last enrichments using the genesets that come out of the first enrichments with the 
biological levels and the ppi's.
  - def __init__(self, gene_list) --> calls all the above functions and puts the results in lists. 
  - def enrich_BL(self) --> does the enrichments by looping through the list of results that is made in the function above. 
- ProteinSet(object) --> contains the code for the actual enrichment. You can give this class a dictionary with terms and 
genesets. 
  - def __init__(self, proteindict, database) --> defining the variables that are going to be used in the whole dataset. 
  - def enrich(self, otherset, background) --> contains the code that loops through the dictionaries to compare the sets of 
  genes. This function calls the function: set_enrichment(). 
  - set_enrichment(self, your_set, other_set, universe, abcd_values=False) --> enriches the genesets. 


Functions for getting the last ranking: 
- FirstRanking() --> Out of each enrichment comes a list of drugs with p-values. These p-values will determine the ranking of the drugs. The smaller the p-value the better (so a drug with a small p-value will be high up in the ranking). Each list will be written to a file. So you will end up with 7 different files. This way, you will be able to take the ranking of each drug and take an average which will be important in the last function (CalLastRanking()).  

Functions for getting the last ranking:  
- LastRanking() --> The list will be ranked by the average ranking.
