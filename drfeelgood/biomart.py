import pandas as pd

class Biomart(object):
    
    def __init__(self):
        """
        Function: 
        ----------
        Making variables that are going to be used throughout the whole code. 
        
        Variables: 
        ----------
        self._biomart_data = a dataset containing ensembl gene id's, protein id's and transcript id's. It also con-
        tains Entrez gene id's, Uniprot id's and the name of the gene. 
        self._idx = a dictionary that converts the dataframe to a readable dictionary. 
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
        k = contains the column that the list starts from and the column that the list should become. 
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
        e = the protein list that is turned into a list of entrez genes. 
        """
        e = self._get_index("protein","Entrez").get(p,None)
        return None if pd.isna(e) else int(e)
    
    def entrez_to_protein(self, e):
        """
        Function: 
        ----------
        Turns a Entrez list of gene id's into a list of ensembl protein id's. 
        """
        return self._get_index("Entrez","protein").get(int(e),None)

