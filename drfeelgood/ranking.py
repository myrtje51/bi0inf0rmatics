from functools import reduce
import pandas as pd

class Ranking(object):
    def __init__(self, enrich_ME):
        """
        Function: 
        ----------
        Making a global variable of the variables that are gonna be used through the whole class. 
        
        Variables: 
        ----------
        self.enrich_ME = de enrichment results from the GiveMeTheDrugs() class, the user can give this as input if he/she
        wants to get back a ranking. 
        """
        self.enrich_ME = enrich_ME
        
    def first_ranking(self):
        """
        Function: 
        ----------
        Ranking all the different biological levels. 
        
        Variables:
        ----------
        only_KEGG = the dataframe (from the enrichment) filtered to be only from the KEGG database.
        only_GO_MF = the dataframe (from the enrichment) filtered to be only from the GO: molecular functions database. 
        only_GO_CC = the dataframe (from the enrichment) filtered to be only from the GO: cellular component database. 
        only_GO_BP = the dataframe (from the enrichment) filtered to be only from the GO: biological process database. 
        only_reactome = the dataframe (from the enrichment) filtered to be only from the reactome database. 
        only_gene_list = the dataframe (from the enrichment) filtered to be only from the given gene list. 
        only_string = the dataframe (from the enrichment) filtered to be only from STRING (PPI's). 
        list_w_ranking = a list with all the ranked enrichment results. 
        """
        only_KEGG = self.enrich_ME[self.enrich_ME['Database'] == "KEGG"]
        only_GO_MF = self.enrich_ME[self.enrich_ME['Database'] == "GO_MF"]
        only_GO_CC = self.enrich_ME[self.enrich_ME['Database'] == "GO_CC"]
        only_GO_BP = self.enrich_ME[self.enrich_ME['Database'] == "GO_BP"]
        only_reactome = self.enrich_ME[self.enrich_ME['Database'] == "Reactome"]
        only_gene_list = self.enrich_ME[self.enrich_ME['Database'] == "GenesOfInterest"]
        only_string = self.enrich_ME[self.enrich_ME['Database'] == "String"]
        
        only_KEGG['Ranking_KEGG'] = only_KEGG['pvalue'].rank()
     
        only_GO_MF['Ranking_GO_MF'] = only_GO_MF['pvalue'].rank()
  
        only_GO_CC['Ranking_GO_CC'] = only_GO_CC['pvalue'].rank()
  
        only_GO_BP['Ranking_GO_BP'] = only_GO_BP['pvalue'].rank()
  
        only_reactome['Ranking_reactome'] = only_reactome['pvalue'].rank()
 
        only_gene_list['Ranking_genes'] = only_gene_list['pvalue'].rank()

        only_string['Ranking_PPI'] = only_string['pvalue'].rank()
     
        
        list_w_ranking = [only_KEGG, only_GO_MF, only_GO_CC, only_GO_BP, only_reactome, only_gene_list, only_string]
        
        return list_w_ranking
    
    def final_ranking(self):
        """
        Function:
        ----------
        Calculating an average ranking from the different rankings made in the function: first_ranking(). 
        
        Variables:
        ----------
        res = the results of the first ranking. 
        for_final = a list with just the names of the drugs in each ranking and the corresponding ranking. 
        name = the drug name.
        df_final = the final dataframe containing all the different rankings and the average ranking of the drugs. This dataframe 
        is sorted by the average ranking. This variable is also returned by the function. 
        """
        res = self.first_ranking()
        
        for_final = []
        for x in res:
            name = x.iloc[:,9:11] 
            for_final.append(name)
            
        df_final = reduce(lambda left,right: pd.merge(left,right,on='drug'), for_final)
        
        df_final['ranking_avg'] = df_final[['Ranking_KEGG', 'Ranking_GO_MF', 'Ranking_GO_CC', 'Ranking_GO_BP', 'Ranking_reactome', 'Ranking_genes', 'Ranking_PPI']].mean(axis=1)

        df_final = df_final.sort_values(by=['ranking_avg'])
    
        return df_final
