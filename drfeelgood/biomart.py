import pandas as pd

class Biomart(object):
    
    def __init__(self):
        """
        Function: 
        ----------
        Making variables that are going to be used throughout the whole code. 
         
        """
        self._biomart_data = pd.read_csv("biomart.tsv", 
              sep='\t', 
              names=["gene", "transcript", "protein", "Entrez", "Uniprot", "name"])
        self._biomart_data = self._biomart_data.applymap(lambda x: None if pd.isna(x) else x)
        
        self._idx = {}
    
    def _get_index(self, col_from, col_to):
        """
        Function:
        ----------
        Checks the index of the columns that are going to be used.  
        
        Variables: 
        ----------
        col_from = from what it needs to be translated.
        col_to   = to what it needs to be translated. 
        
        Results: 
        ----------
        self._idx = dictionary with the translation. 
        """
        k = (col_from, col_to)
        if k not in self._idx:
            self._idx[k] = { r[col_from] : r[col_to] for i,r in self._biomart_data.iterrows() }
        
        return self._idx[k]

    def protein_to_entrez(self, p):
        """
        Function: 
        ----------
        Turns a list of ensembl protein id's into a list of Entrez gene id's. 
        
        Variables: 
        ----------
        p = a protein that needs to be translated. 
        
        Returns:
        ---------
        None is it can't find a translation, else the translation (entrez id). 
        """
        e = self._get_index("protein","Entrez").get(p,None)
        return None if pd.isna(e) else int(e)
    
    def entrez_to_protein(self, e):
        """
        Function: 
        ----------
        Turns a Entrez list of gene id's into a list of ensembl protein id's. 
        
        Variables: 
        ----------
        e = an entrez id that needs to be translated. 
        
        Returns:
        ----------
        the translation(s) of the entrez id to protein id's (ensembl). 
        """
        return self._get_index("Entrez","protein").get(int(e),None)
    
    def name_to_entrez(self, n):
        """
        Function: 
        -----------
        Turns a gene name into an entrez gene id. 
        
        Variabeles: 
        -----------
        n = a gene name that needs to be translated.
        
        Returns:
        -----------
        None if there is no translation, else the translation (entrez id). 
        """
        e = self._get_index("name", "Entrez").get(n,None) 
        return None if pd.isna(e) else int(e) 

