import numpy as np
import pandas as pd

class RocAuc(object): 
    
    def readfilter_databases(self, source):
        """
        Function:
        ----------
        This function reads in the DrugAge dataset. This is a dataset that contains a bunch of drugs that are known to increase or
        decrease a lifespan. It also reads in the mapped dataset with drugs and targets. It then filters the DrugAge dataset by 
        only keeping the drugs that are also in the mapped dataset. 
        
        Variables:
        ----------
        source = the database that was given by the user. 
        
        Returns:
        ----------
        source = the database that was given by the user. 
        
        """
        mapped = pd.read_csv('~/drfeelgood/Files/mapped_DB_STITCH_actions_first.tsv', sep='\t')
   
        DT_list = set(mapped['Name'].unique())
    
        DADT_pro = source & DT_list
        
        return source
    
    def ROC(self, df, rank):
        """
        Function:
        ----------
        This function returns a list that can be used in matplotlib to make a ROC plot. 
        
        Variables:
        ----------
        df      = a dataframe with a ranking.
        rank    = a name of the database that was used to make the ranking. 
        
        Returns: 
        ---------
        Two lists of ROC's. 
        """
        df = df.copy()
        df.ranking = df[rank]
        roc = []
        for t in sorted(np.unique(df.ranking.values)):
            tp = df[(df.ranking <= t) & df.TP].shape[0]
            tn = df[(df.ranking <= t) & ~df.TP].shape[0]
            fp = df[(df.ranking > t) & ~df.TP].shape[0]
            fn = df[(df.ranking > t) & df.TP].shape[0]
            
            tpr = tp / (tp + fn)
            fpr = fp / (fp + tn)
            roc.append((1-fpr,tpr))
        
        new_roc = [ roc[0] ]
        for fpr, tpr in roc[1:]:
            if new_roc[-1][0] == fpr:
                new_roc[-1] = (fpr, max(tpr, new_roc[-1][1]))
            else:
                new_roc.append((fpr, tpr))
        
        return list(zip(*new_roc))
    
    def AUC(self, fpr, tpr, topn=1):
        """
        Function: 
        ----------
        This function calculates the AUC (Area Under the Curve) 
        
        Variables: 
        ----------
        fpr = False Positive Rate.
        tpr = True Positive Rate. 
        
        Returns: 
        ---------
        pauc = The (partial) Area Under the Curve. 
        
        """
        auc = 0.0
        newtpr = 1.0
        
        fpr = [0] + list(fpr)
        tpr = [0] + list(tpr)
        
        for i in range(len(fpr)-1):
            i = i+1
            if fpr[i] >= topn:
                newtpr = tpr[i]
                auc += (topn-fpr[i-1]) * tpr[i-1]
                break
            auc += (fpr[i] - fpr[i-1]) * tpr[i-1]
    
        surface = newtpr*topn
       
        if surface > 0:
            pauc = auc/float(surface) 
        else: 
            pauc = 0
    
        return pauc
