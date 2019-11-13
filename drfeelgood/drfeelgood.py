from .proteinset import ProteinSet
from .biomart import Biomart
bm = Biomart() 
from .databases import Databases
db = Databases() 
from rpy2.robjects.vectors import StrVector
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects import r, pandas2ri
import pandas as pd
import numpy as np
import re


class GiveMeTheDrugs(object):
    def __init__(self, gene_list):
        """
        Function: 
        ----------
        This function calls all the other functions that have something to do with enriching with R. It also imports the needed
        R-packages. Before this function can be called please make sure you have the packages already installed and the packages 
        are in the right folder (/Users/user/Anaconda3/Lib/R/library). 

        Variables: 
        ----------
        self.ensembl = all the target proteins from the STITCH database. 
        self.get_entr_filtered_ens = all the ensembl proteins turned into gene ID's.
        acr_res = dictionary of gene_list: and a list of genes. 
        gene_set = a list of ageing related genes.
        clusterProfiler = the package: clusterProfiler imported using rpy2.
        ReactomePA = the package: ReactomePA imported using rpy2.
        KEGG_res = the enrichment called where the Drug-Target dataset gets enriched against the first enrichent done. 
        GO_MF_en_res = the enrichment called where the Drug-Target dataset gets enriched against the first enrichment done.
        GO_CC_en_res = the enrichment called where the Drug-Target dataset gets enriched against the first enrichment done.  
        GO_BP_en_res = the enrichment called where the Drug-Target dataset gets enriched against the first enrichment done.  
        reactome_en_res = the enrichment called where the Drug-Target dataset gets enriched against the first enrichment done. 
        PPI_rest = the enrichment called where the Drug-Target dataset gets enriched against the first enrichment done.
        GOI = the enrichment called where the Drug-Target dataset gets enriched against the ageing related genelist.
        self.PS = list with all the enrichments classes called. 
        self.dictio_drugs = calls the function: make_dictio_DT() 
        """
        
        ensembl = pd.read_csv('/home/mhaan/STITCH_proteins.txt')
        ensembl = list(ensembl['protein'])
        self.ensembl = ensembl
        
        get_entr_filtered_ens = list(map(bm.protein_to_entrez, self.ensembl))
        self.get_entr_filtered_ens = get_entr_filtered_ens
        self.get_entr_filtered_ens = [i for i in self.get_entr_filtered_ens if i]
        
        gene_set1 = list(gene_list)
        gene_set = list(map(int, gene_set1))
        
        get_ensp_filtered = list(map(bm.entrez_to_protein, gene_set))
        get_ensp_filtered = [i for i in get_ensp_filtered if i]

        clusterProfiler = rpackages.importr('clusterProfiler')
        ReactomePA = rpackages.importr('ReactomePA')
        GO_database = rpackages.importr('org.Hs.eg.db') 
        
        print("prep kegg")
        KEGG_res = db.cluster_profiler_KEGG(clusterProfiler, gene_list)
        print("prep go_mf")
        GO_MF_en_res = db.cluster_profiler_GO_MF(clusterProfiler, gene_list, GO_database)
        print("prep go_cc")
        GO_CC_en_res = db.cluster_profiler_GO_CC(clusterProfiler, gene_list, GO_database)
        print("prep go_bp")
        GO_BP_en_res = db.cluster_profiler_GO_BP(clusterProfiler, gene_list, GO_database)
        print("prep reactome")
        reactome_en_res = db.Reactome(clusterProfiler, ReactomePA, gene_list)
        print("prep PPI's")
        ppi_res = db.make_dictio_ppi(get_ensp_filtered, self.get_entr_filtered_ens)
        print("prep genes of interest")
        GOI = ProteinSet({ "gene_list" : list(map(int, gene_list)) }, "GenesOfInterest", set(self.get_entr_filtered_ens))
        self.PS = [KEGG_res, GO_MF_en_res, GO_CC_en_res, GO_BP_en_res, reactome_en_res, GOI, ppi_res]

        self.dictio_drugs = db.make_dictio_DT()

    
    def enrich_BL(self):
        """
        Function: 
        ----------
        This function does all the enrichments that were previously called with the function ProteinSet(object).
        
        Variables:
        ----------
        super_x = list with the top results for every drug in every database. 
        df = the dataframe with the enrichment results (the really low results are in there as well)
        self.drug_enrichments = all the top results for every drug in every database in a dataframe. 
        """

        super_x = []
        for i, (drugs, targets) in enumerate(self.dictio_drugs.items()):
            print("%d/%d: %s" % (i+1, len(self.dictio_drugs), drugs))
            df = pd.concat([ps.enrich(targets, set(self.get_entr_filtered_ens)) for ps in self.PS])
            if len(df) != 0:
                df['drug'] = drugs
                df = df.sort_values("pvalue").groupby("Database", as_index=False).first() 
            super_x.append(df)
        self.drug_enrichments = pd.concat(super_x)
        return self.drug_enrichments
    
