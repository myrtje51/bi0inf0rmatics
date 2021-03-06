from collections import namedtuple
import scipy.stats as sstats
import statsmodels.stats.multitest as multi
import pandas as pd

class ProteinSet(object):
    def __init__(self, proteindict, database, dataset_universe):
        """
        Function: 
        ----------
        Making a global variable of the variables that are gonna be used through the whole class. 
        
        Variables: 
        ----------
        self.proteindict = a dictionary with a term and a list of proteins per item. 
        self.database    = the name of the database. 
        """
        self.proteindict = { name : set(p) for name, p in proteindict.items() }
        self.database = database
        self.dataset_universe = dataset_universe
    
    def enrich(self, otherset, background):
        """
        Function: 
        ----------
        This function takes 3 sets of proteins (or genes) and uses them to make an enrichment using either fisher's exact test or 
        the chi2 (depending on how big the sets are). 
        
        Variables:
        ----------
        list_res  = a list with lists that will later be turned into a dataframe. Eech list within the list will have information
                    about a row in the table. 
        name      = the name of the drug 
        pset      = a set of proteins that are targets of the drug. 
        term      = the name in a list. 
        proteins  = the proteins in a list. 
        results   = the enrichment results in a NamedTuple. 
        l_results = the enrichment results turned into a list. 
        joined    = the name of the database and the term merged with the l_results list. 
        df_final  = a dataframe with all the enrichment results. 
        """
        list_res = []
        
        dataset_universe_background = background & self.dataset_universe
        
        for name, pset in self.proteindict.items():
            if len(set(pset) & set(otherset)) == 0:
                continue
            results = self.set_enrichment(pset, otherset, dataset_universe_background)
            joined = [self.database, name] + list(results) + [ list(pset) ]
            list_res.append(joined)
        
        if len(list_res) != 0: 
            df_final = pd.DataFrame(list_res)
            df_final.columns = ['Database', 'Name', 'oddsratio', 'c2statistic', 'pvalue', 'table', 'method', 'proteins'] 
            padjust = multi.multipletests(df_final['pvalue'], method='fdr_by')
            padjust = padjust[1][0]
            df_final['padjust'] = padjust
        else: 
            df_final = pd.DataFrame(list_res)
        return df_final
    
    def set_enrichment(self, your_set, other_set, universe, abcd_values=False):
    
        resTuple = namedtuple("setEnrichmentResult", [ 'oddsratio', 'c2statistic', 'pvalue', 'table', 'method'])

        universe  = set(universe)
        your_set  = set(your_set) & universe
        other_set = set(other_set) & universe

        a = your_set & other_set
        b = other_set - your_set
        c = your_set - other_set
        d = universe - (your_set | other_set)

        table = [ [len(a), len(b)], [len(c), len(d)]]
        
        method = 'fisher'
        oddsratio, p = sstats.fisher_exact(table)
        chi2 = None

        if abcd_values:
            return resTuple(oddsratio, chi2, p, [[a,b],[c,d]], method)
        else:
            return resTuple(oddsratio, chi2, p, table, method)


