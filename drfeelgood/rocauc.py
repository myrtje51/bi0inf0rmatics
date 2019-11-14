import numpy as np
import pandas as pd

class RocAuc(object): 
    
    def readfilter_databases(self):
        """
        Function:
        ----------
        This function reads in the DrugAge dataset. This is a dataset that contains a bunch of drugs that are known to increase or
        decrease a lifespan. It also reads in the mapped dataset with drugs and targets. It then filters the DrugAge dataset by 
        only keeping the drugs that are also in the mapped dataset. 
        
        Variables:
        ----------
        drugAge = the DrugAge dataset. 
        mapped = the mapped drug target dataset. 
        pro_longevity = all the pro-longevity drugs in the DrugAge dataset. 
        anti_longevity = all the anti_longevity drugs in the DrugAge dataset.
        DA_list_anti = anti-longevity but made into a set.
        DA_list_pro = pro-longevity but made into a set. 
        DADT_pro = all the pro-longevity drugs present in both the DrugAge dataset and the mapped dataset.
        
        """
        drugAge = pd.read_csv('drugage.csv')
        mapped = pd.read_csv('mapped_DB_STITCH_actions_first.tsv', sep='\t')
    
        pro_longevity = drugAge[drugAge['avg_lifespan_change'] > 0]
        anti_longevity = drugAge[drugAge['avg_lifespan_change'] < 0]
        
        DA_list_anti = set(anti_longevity['compound_name'].unique())
        DA_list_pro = set(pro_longevity['compound_name'].unique())
        DT_list = set(mapped['Name'].unique())
    
        DADT_pro = DA_list_pro & DT_list
        
        return DA_list_pro
    
    def ROC(self, df, rank):
        """
        Function:
        ----------
        This function returns a list that can be used in matplotlib to make a ROC plot. 
        
        Variables:
        ----------
        df = a dataframe with a ranking.
        roc = a list with points for a plot. 
        tp = true postives
        tn = true negatives
        fp = false positives
        fn = false negatives
        tpr = true positive rate
        fpr = false positive rate 
        new_roc = to avoid two of the same tpr/fpr right after the other. 
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
    
    def AUC(self, fpr, tpr):
        """
        Function: 
        ----------
        This function calculates the AUC (Area Under the Curve) 
        
        Variables: 
        ----------
        auc = the Area Under the Curve. 
        
        """
        auc = 0
        for i in range(len(fpr)-1):
            auc += (fpr[i+1] - fpr[i]) * tpr[i]
        return auc