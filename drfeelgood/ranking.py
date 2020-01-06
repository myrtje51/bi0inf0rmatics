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
        enrich_ME = de enrichment results from the GiveMeTheDrugs() class, the user can give this as input if he/she
                         wants to get back a ranking. 
        """
        self.enrich_ME = enrich_ME
        
    def ranking1(self):
        """
        Function: 
        ----------
        Ranking all the different biological levels. 
        
        Returns:
        ----------
        list_w_ranking = list with the different rankings. 
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
        
        Returns:
        ----------
        df_final = the final dataframe with the ranked drugs for every biological level. 
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

