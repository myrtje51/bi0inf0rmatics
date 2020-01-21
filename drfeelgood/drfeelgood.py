from .proteinset import ProteinSet
from .biomart import Biomart
bm = Biomart() 
from .ranking import Ranking
from .rocauc import RocAuc
from .databases import Databases
db = Databases() 
from rpy2.robjects.vectors import StrVector
import rpy2.robjects.packages as rpackages
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robjects
import matplotlib.pylab as plt
import pandas as pd
import numpy as np
import re
from fuzzywuzzy import fuzz
from fuzzywuzzy import process


class GiveMeTheDrugs(object):
    def __init__(self, gene_list, db_list=None, db_name=None):
        """
        Function: 
        ----------
        This function calls all the other functions that have something to do with enriching with R. It also imports the needed
        R-packages. Before this function can be called please make sure you have the packages already installed and the packages 
        are in the right folder (/Users/user/Anaconda3/Lib/R/library). 

        Variables: 
        ----------
        gene_list = list of genes the user wants to find drugs for. 
        db_list   = list with biological levels that the user wants to use. 
        db_name   = Name of the drug-target database that the user wants to use (CTD or DrugBank). 
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
        
        if db_list is None: 
            db_list = ['bp', 'cc', 'mf', 'Genes', 'KEGG', 'reactome', 'string']
        
        self.database_list = [x.lower() for x in db_list]
        
        self.PS = []
        
        if "kegg" in self.database_list: 
            print("prep kegg")
            KEGG_res = db.cluster_profiler_KEGG(clusterProfiler, gene_list)
            self.PS.append(KEGG_res)
        if "go_mf" in self.database_list or "molecular function" in self.database_list or "mf" in self.database_list: 
            print("prep go_mf")
            GO_MF_en_res = db.cluster_profiler_GO_MF(clusterProfiler, gene_list, GO_database)
            self.PS.append(GO_MF_en_res)
        if "go_cc" in self.database_list or "cellular component" in self.database_list or "cc" in self.database_list: 
            print("prep go_cc")
            GO_CC_en_res = db.cluster_profiler_GO_CC(clusterProfiler, gene_list, GO_database)
            self.PS.append(GO_CC_en_res)
        if "go_bp" in self.database_list or "biological process" in self.database_list or "bp" in self.database_list: 
            print("prep go_bp")
            GO_BP_en_res = db.cluster_profiler_GO_BP(clusterProfiler, gene_list, GO_database)
            self.PS.append(GO_BP_en_res)
        if "reactome" in self.database_list:
            print("prep reactome")
            reactome_en_res = db.Reactome(clusterProfiler, ReactomePA, gene_list)
            self.PS.append(reactome_en_res)
        if "string" in self.database_list: 
            print("prep PPI's")
            ppi_res = db.make_dictio_ppi(get_ensp_filtered, self.get_entr_filtered_ens)
            self.PS.append(ppi_res)
        if "genes" in self.database_list: 
            print("prep genes of interest")
            GOI = ProteinSet({ "gene_list" : list(map(int, gene_list)) }, "GenesOfInterest", set(self.get_entr_filtered_ens))
            self.PS.append(GOI)

        if db_name is None: 
            print("User didn't specify which database to use. Now using the default database: DrugBank...")
            self.dictio_drugs = db.make_dictio_DT()
        elif db_name == "CTD" or db_name == "Comparative Toxicogenomics Database": 
            self.dictio_drugs = db.get_CTD() 
        elif db_name == "DB" or db_name == "DrugBank" or db_name == "drugbank":
            self.dictio_drugs = db.make_dictio_DT() 

    
    def DrugDb_enrich(self):
        """
        Function: 
        ----------
        This function does all the enrichments that were previously called with the function ProteinSet(object).
        
        Returns:
        ----------
        self.drug_enrichments = a dataframe with all the enrichments that were done (unranked). 
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
    
    def first_ranking(self, enrichment):
        """
        Function:
        ----------
        This function gives each drug a ranking from each database. This is still without the average ranking. 
        
        Variables:
        ----------
        enrichment = a dataframe with the enrichment results.
        
        Returns: 
        ----------
        The results of the function: ranking1().
        """
        r = Ranking(enrichment)
        return r.ranking1() 
    
    def final_ranking(self, enrichment):
        """
        Function:
        ----------
        This function calculates the average ranking of each drug based on the enrichment that is done on each biological level.
        
        Variables:
        ----------
        enrichment = a dataframe with the enrichment results. 
        
        Returns:
        ----------
        self.FR = the results of ranking2() a.k.a. the final ranking. 
        """
        r = Ranking(enrichment)
        self.FR = r.ranking2() 
        return self.FR
    
    def rocauc_maker(self, source, topn=None, FiR=None):
        """
        Function:
        ----------
        This function calculates the AUC and plots a ROC. 
        
        Variables: 
        ----------
        source   = list of drugs that the program has to compare the predictions with. 
        topn     = an integer if the user wants to only use a portion of the prediction (like the top 100 predicted drugs) 
        FiR      = the final ranking. 
        
        """
        RA = RocAuc() 
        
        if isinstance(FiR, type(None)):
            FiR = self.FR
            
        names = list(FiR.columns)
        
        if 'TP' in names:
            names.remove('TP')
         
        true_pos = RA.readfilter_databases(source)
        TrP = []
        for d in list(FiR.drug):
            fuzzystring = process.extract(d, true_pos) 
            best = [item for item in fuzzystring if item[1] == 100]
            if len(best) > 0: 
                TrP.append(True)
            elif len(best) == 0: 
                TrP.append(False) 
        
        FiR['TP'] = TrP
        
        if topn is None:
            topn = 1
        
        roc = [ RA.ROC(FiR.sort_values(r), r) for r in names[1:] ]
        auc = { r: RA.AUC(roc[i][0], roc[i][1], topn) for i,r in enumerate(names[1:]) }
        
        fig = plt.figure(figsize=(10,10), dpi=300)
        ax = plt.gca()
        for i,r in enumerate(names[1:]):
            ax.plot([ name for name in roc[i][0] if name <= topn ], roc[i][1][:len([ name for name in roc[i][0] if name <= topn ])], label='%s: %f' % (r, auc[r]))
        x_lims = ax.get_xlim()
        y_lims = ax.get_ylim()
        ax.plot(x_lims,y_lims, c='r')
        ax.legend()
        return auc
