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
        
    def ranking1(self):
        """
        Function: 
        ----------
        Ranking all the different biological levels. 
        
        Variables:
        ----------
        only_KEGG      = the dataframe (from the enrichment) filtered to be only from the KEGG database.
        only_GO_MF     = the dataframe (from the enrichment) filtered to be only from the GO: molecular functions database. 
        only_GO_CC     = the dataframe (from the enrichment) filtered to be only from the GO: cellular component database. 
        only_GO_BP     = the dataframe (from the enrichment) filtered to be only from the GO: biological process database. 
        only_reactome  = the dataframe (from the enrichment) filtered to be only from the reactome database. 
        only_gene_list = the dataframe (from the enrichment) filtered to be only from the given gene list. 
        only_string    = the dataframe (from the enrichment) filtered to be only from STRING (PPI's). 
        list_w_ranking = a list with all the ranked enrichment results. 
        """
        list_w_ranking = []
        for d in self.enrich_ME.Database.unique():
            only = self.enrich_ME[self.enrich_ME['Database'] == d]
            only[d] = only['pvalue'].rank()
            list_w_ranking.append(only) 
        
        return list_w_ranking
    
    def ranking2(self):
        """
        Function:
        ----------
        Calculating an average ranking from the different rankings made in the function: first_ranking(). 
        
        Variables:
        ----------
        res       = the results of the first ranking. 
        for_final = a list with just the names of the drugs in each ranking and the corresponding ranking. 
        name      = the drug name.
        df_final  = the final dataframe containing all the different rankings and the average ranking of the drugs. This dataframe 
                    is sorted by the average ranking. This variable is also returned by the function. 
        """
        res = self.ranking1()
        
        for_final = []
        for x in res:
            name = x.iloc[:,9:11] 
            for_final.append(name)
            
        df_final = reduce(lambda left,right: pd.merge(left,right,on='drug'), for_final)
        
        df_final['ranking_avg'] = df_final[[d for d in self.enrich_ME.Database.unique()]].mean(axis=1) 
     
        df_final = df_final.sort_values(by=['ranking_avg'])
    
        return df_final
