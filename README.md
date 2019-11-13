# drfeelgood

### For new users, some packages that are needed to run this program:

*Python packages:* 
- [rpy2](https://anaconda.org/r/rpy2)
- [pandas](https://anaconda.org/anaconda/pandas)
- [numpy](https://anaconda.org/anaconda/numpy)
- [scipy](https://anaconda.org/anaconda/scipy)
- [requests](https://anaconda.org/anaconda/requests)
- [statsmodels](https://pypi.org/project/statsmodels/)

*R packages (within anaconda/miniconda):*
- [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
- [ReactomePA](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html) 
- [STRINGdb](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)
- [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

### What does the package do (summary)? 
Drfeelgood is a package that takes a list of genes (preferably an R-array which can be read from a file) and enriches this list of genes against a few different biological levels (KEGG, reactome, GO and protein-protein interactions). It then enriches the lists that come out of the first enrichments against a dataset containing drugs with their corresponding targets. This enrichment is done with only the fisher's exact test. 
These enrichments are then ranked by their pvalue. When the lists have been ranked an average ranking is calculated. A list of drugs will be returned sorted by the average ranking. 

