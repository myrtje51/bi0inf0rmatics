import pandas as pd
import numpy as np
import scipy as sp
from drfeelgood import biomart 
bm = biomart.Biomart()
import scipy.linalg
import networkx as nx
import matplotlib.pyplot as plt 


class Diffusion(object):
    def __init__(self, genes):
        """
        Variables: 
        ------------
        genes = A list of genes that the user wants to invesigate. 
        
        """
        self.genes = genes
        
        string = pd.read_csv("9606.protein.links.v11.0.txt", sep=' ')
        string = string[string.combined_score >= 700]
        
        self.gn = list(set(string.protein1) | set(string.protein2))
        
        dt = pd.read_csv("mapped_DB_STITCH_actions_first.tsv", sep='\t')
        dt = dt[dt['item_id_b'].isin(list(string.protein1)) | dt['item_id_b'].isin(list(string.protein2))]
        
        self.dn = list(set(dt.Name))
        
        ndrugs = len(set(dt.Name))
        ngenes = len(set(string.protein1))
        self.matrix = np.zeros((ngenes+ndrugs+1, ngenes+ndrugs+1), dtype=float)
        
    def disease(self):
        """
        Returns:
        ----------
        self.matrix = the matrix with protein-protein interactions, drugs and the genes that are related to 
                      the disease the user wants to investigate. 
        """

        GM = {n:i for i, n in enumerate(self.gn)}
        
        nb_drugs = len(GM)
        end = len(self.dn) + nb_drugs
        l_dg = []
        for x in range(len(GM), end):
            l_dg.append(x)
        DM = {i:n for i, n in zip(self.dn, l_dg)}
        
        protein_ids = list(set(map(bm.entrez_to_protein, self.genes)))
        protein_ids = [i for i in protein_ids if i]
        protein_h = ['9606.' + p for p in protein_ids]
        protein_h = set(protein_h) & set(self.gn)
        
        self.matrix = self.fill(DM, GM, protein_h)
        
        return self.matrix
    
    def fill(self, DM, GM, protein_h):
        """
        Variables:
        ------------
        DM = A dictionary with drugs and the corresponding index in the matrix. 
        GM = A disctionary with genes and the corresponding index in the matrix. 
        protein_h = a set of proteins related to the disease that needs to be investigated. 
        
        Returns: 
        ------------
        self.matrix = the matrix with protein-protein interactions, drugs and the genes that are related to 
                      the disease the user wants to investigate. 
        """
        
        for i,r in dt.iterrows():
            self.matrix[DM[r.Name], GM[r.item_id_b]] = 1
            self.matrix[GM[r.item_id_b], DM[r.Name]] = 1
        
        for i,r in string.iterrows():
            norm = r.combined_score/1000
            self.matrix[GM[r.protein1], GM[r.protein2]] = norm
            self.matrix[GM[r.protein2], GM[r.protein1]] = norm
        
        for r in protein_h:
            self.matrix[-1, GM[r]] = 1
            self.matrix[GM[r], -1] = 1
        
        return self.matrix
    
    def zeros(self):
        """
        Returns:
        -----------
        S = a numpy array filled with zeros, where the disease gets a value of -1. 
        """
        S = np.zeros(self.matrix.shape[0])
        S[-1] = 1
        return S
    
    def laplacian(self, beta=0.01):
        """
        Variables:
        ------------
        beta = how many times the matrix needs to be multiplied. (default = 0.01) 
        
        Returns: 
        ------------
        The laplacian matrix is returned. 
        """
        return sp.linalg.expm(beta*scipy.sparse.csgraph.laplacian(self.matrix))
    
    def diffusion(self):
        """
        Returns:
        ----------
        P = The diffusion results with the heatflow. 
        """
        S = self.zeros()
        D = self.laplacian(beta=0.002)
        P = np.dot(S, D)
        return P
        
    def ranking(self, P):
        """
        Variables:
        -----------
        P = The diffusion results with the heatflow. 
        
        Results:
        -----------
        ranked = A ranked dataframe with the drugs ranked based on the heatflow. 
        """
        nodes = self.gn + self.dn + ["Disease"]
        heat_flow = list(zip(nodes, P))
        
        df = pd.DataFrame(heat_flow)
        df.columns = ['Nodes', 'HeatFlow']
        
        drugsdf = df[df['Nodes'].isin(self.dn)]
        
        drugsdf['Ranking'] = drugsdf['HeatFlow'].rank(ascending=False)
        
        ranked = drugsdf.sort_values(by=['Ranking'], ascending=True)
        
        return ranked     
