from .proteinset import ProteinSet
from .biomart import Biomart 
bm = Biomart() 
from rpy2.robjects.vectors import StrVector
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects import r, pandas2ri
import pandas as pd
import re

class Databases(object):
    
    def make_dictio_DT(self): 
        """
        Function: 
        ----------
        This function makes a dictionary with a drug and a list of proteins that are targets of that drug.
        
        Returns: 
        ----------
        dictio = a dictionary with a drug and the corresponding targets (proteins). Where the drug has two or more targets. 
        """
        mapped = pd.read_csv('mapped_DB_STITCH_actions_first.tsv', sep='\t')
        mapped['item_id_b'] = mapped['item_id_b'].map(lambda x: x.lstrip('9606.'))
        mapped = mapped[['CID', 'InChIKey', 'DrugBank ID', 'Name', 'item_id_b']].drop_duplicates()
        
        mapped['Entrez'] = list(map(bm.protein_to_entrez, list(mapped['item_id_b'])))
        mapped = mapped.dropna(subset=['Entrez'])
        mapped['Entrez'] = mapped['Entrez'].astype(int)
        
        dictio = {}
        for i in mapped['Name'].unique(): 
            dictio[i] = [mapped['Entrez'][j] for j in mapped[mapped['Name']==i].index]
        
        dictio = {k: v for k, v in dictio.items() if len(v) > 2}
        
        return dictio 
    
    def get_CTD(self):
        """
        Function:
        ----------
        This function turns the CTD database into a dictionary with drugs and their gene targets (names)
        
        Returns: 
        ----------
        drug_genes = a dictionary with a drug and the corresponding targets (proteins). Where the drug has two or more targets.
        
        """
        chemgene = pd.read_csv('CTD_chem_gene_ixns.tsv', sep='\t', comment='#', names=['ChemicalName','ChemicalID','CasRN','GeneSymbol','GeneID','GeneForms','Organism','OrganismID','Interaction','InteractionActions','PubMedIDs'])
        chemgeneH = chemgene[chemgene['Organism'] == 'Homo sapiens']
        chemgeneH = chemgeneH.drop(columns=['ChemicalID', 'CasRN', 'OrganismID', 'PubMedIDs'])
        
        chemdis = pd.read_csv('CTD_chemicals_diseases.tsv', sep='\t', comment='#', names=['ChemicalName','ChemicalID','CasRN','DiseaseName','DiseaseID','DirectEvidence','InferenceGeneSymbol','InferenceScore','OmimIDs','PubMedIDs'])
        chemdisT = chemdis[chemdis['DirectEvidence'] == 'therapeutic']
        chemdisT = chemdisT.drop(columns=['OmimIDs', 'PubMedIDs', 'DiseaseID', 'InferenceGeneSymbol', 'InferenceScore', 'DiseaseName'])
        
        THmerge = pd.merge(chemgeneH, chemdisT, on='ChemicalName')
        THmerge = THmerge.drop_duplicates() 
        
        drug_genes = {}
        for i in THmerge['ChemicalName'].unique(): 
            drug_genes[i] = set([THmerge['GeneID'][j] for j in THmerge[THmerge['ChemicalName']==i].index]) 
        
        drug_genes = {k: v for k, v in drug_genes.items() if len(v) > 2}
        
        return drug_genes  
        
    def cluster_profiler_KEGG(self, clusterProfiler, data): 
        """
        Function: 
        ----------
        This function uses the R-package: clusterProfiler to enrich a list of genes for KEGG pathways. You will get back a list
        of KEGG pathways and the sets of genes part of that pathway. The list is sorted by the adjusted p-value.
        
        Variables: 
        ----------
        clusterProfiler = The R package: clusterProfiler.
        data            = List of genes. 
        
        Returns: 
        ----------
        the call of the class: ProteinSet(). 
        """
        
        enrich_KEGG = clusterProfiler.enrichKEGG(data, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = 'BY')
        KEGGdat = enrich_KEGG.slots['result']
        
        KEGG = clusterProfiler.download_KEGG('hsa', keggType = "KEGG", keyType = "kegg")
        full_gl_KEGG = set(KEGG[0][1])
        full_gl_KEGG = set(map(int, full_gl_KEGG))
        
        df = pd.DataFrame(index=range(len(KEGGdat[0])))
        df['ID'] = KEGGdat[0]
        df['Description'] = KEGGdat[1]
        df['GeneRatio'] = KEGGdat[2]
        df['pvalue'] = KEGGdat[4]
        df['padjust'] = KEGGdat[5]
        df['Count'] = KEGGdat[8]
        empty_lijst_KEGG = []
        for x in KEGGdat[7]:
            per_pathway = []
            lijst = x.split('/')
            lijst = list(map(int, lijst))
            empty_lijst_KEGG.append(lijst)
        df['GeneID'] = empty_lijst_KEGG
        
        dictionary = {} 
        for index, row in df.iterrows():
            dictionary[row['ID']] = row['GeneID']
            
        return ProteinSet(dictionary,"KEGG",full_gl_KEGG)
    
    def cluster_profiler_GO_MF(self, clusterProfiler, data, GO_database): 
        """
        Function: 
        ----------
        This function uses the R-package: clusterProfiler to enrich a list of genes for GO annotations (specifically molecular
        function). You will get back a list of GO annotations with the sets of genes that have something to do with that GO anno-
        tation. This list is sorted by the adjusted p-value. 
        
        Variables: 
        ----------
        clusterProfiler = The R package: clusterProfiler.
        data            = List of genes. 
        GO_database     = database with GO-terms. 
        
        Results: 
        ----------
        the call of the class: ProteinSet(). 
        """
        
        enrich_GO_MF = clusterProfiler.enrichGO(data, 'org.Hs.eg.db', ont = 'MF', pvalueCutoff = 0.05, pAdjustMethod = 'BY')
        MFdat = enrich_GO_MF.slots['result']
        
        robjects.r('k <- keys(org.Hs.eg.db, "ENTREZID")')
        robjects.r('df <- select(org.Hs.eg.db, keys=k, columns=c("ONTOLOGY","GO"),keytype="ENTREZID")')
        MF = robjects.r('df[which(df$ONTOLOGY == "MF"), ]')
        full_gl_GOMF = set(MF[0])
        full_gl_GOMF = set(map(int, full_gl_GOMF))
        
        df_GO_MF = pd.DataFrame(index=range(len(MFdat[0])))
        df_GO_MF['ID'] = MFdat[0]
        df_GO_MF['Description'] = MFdat[1]
        df_GO_MF['GeneRatio'] = MFdat[2]
        df_GO_MF['pvalue'] = MFdat[4]
        df_GO_MF['padjust'] = MFdat[5]
        df_GO_MF['Count'] = MFdat[8]
        empty_lijst_GO_MF = []
        for x in MFdat[7]: 
            lijst = x.split('/')
            lijst = list(map(int, lijst))
            empty_lijst_GO_MF.append(lijst)
        df_GO_MF['GeneID'] = empty_lijst_GO_MF
        
        dictionary = {} 
        for index, row in df_GO_MF.iterrows():
            dictionary[row['ID']] = row['GeneID']
            
        return ProteinSet(dictionary,"GO_MF",full_gl_GOMF)
    
    def cluster_profiler_GO_CC(self, clusterProfiler, data, GO_database): 
        """
        Function: 
        ----------
        This function uses the R-package: clusterProfiler to enrich a list of genes for GO annotations (specifically cellular com-
        ponent). You will get back a list of GO annotations with the sets of genes that have something to do with that GO anno-
        tation. This list is sorted by the adjusted p-value. 
        
        Variables: 
        ----------
        clusterProfiler = The R package: clusterProfiler.
        data            = List of genes. 
        GO_database     = database with GO-terms. 
        
        Results:
        ----------
        the call of the class: ProteinSet(). 
        """
        
        enrich_GO_CC = clusterProfiler.enrichGO(data, 'org.Hs.eg.db', ont = 'CC', pvalueCutoff = 0.05, pAdjustMethod = 'BY')
        CCdat = enrich_GO_CC.slots['result']
        
        robjects.r('k <- keys(org.Hs.eg.db, "ENTREZID")')
        robjects.r('df <- select(org.Hs.eg.db, keys=k, columns=c("ONTOLOGY","GO"),keytype="ENTREZID")')
        CC = robjects.r('df[which(df$ONTOLOGY == "CC"), ]')
        full_gl_GOCC = set(CC[0])
        full_gl_GOCC = set(map(int, full_gl_GOCC))
        
        df_GO_CC = pd.DataFrame(index=range(len(CCdat[0])))
        df_GO_CC['ID'] = CCdat[0]
        df_GO_CC['Description'] = CCdat[1]
        df_GO_CC['GeneRatio'] = CCdat[2]
        df_GO_CC['pvalue'] = CCdat[4]
        df_GO_CC['padjust'] = CCdat[5]
        df_GO_CC['Count'] = CCdat[8]
        empty_lijst_GO_CC = []
        for x in CCdat[7]: 
            lijst = x.split('/')
            lijst = list(map(int, lijst))
            empty_lijst_GO_CC.append(lijst)
        df_GO_CC['GeneID'] = empty_lijst_GO_CC
        
        dictionary = {} 
        for index, row in df_GO_CC.iterrows():
            dictionary[row['ID']] = row['GeneID']
            
        return ProteinSet(dictionary,"GO_CC",full_gl_GOCC)
    
    def cluster_profiler_GO_BP(self, clusterProfiler, data, GO_database):
        """
        Function: 
        ----------
        This function uses the R-package: clusterProfiler to enrich a list of genes for GO annotations (specifically biological
        process). You will get back a list of GO annotations with the sets of genes that have something to do with that GO anno-
        tation. This list is sorted by the adjusted p-value.
        
        Variables: 
        ----------
        clusterProfiler = The R package: clusterProfiler.
        data            = List of genes. 
        GO_database     = database with GO-terms. 
        
        Results: 
        ----------
        the call of the class: ProteinSet(). 
        """
        
        enrich_GO_BP = clusterProfiler.enrichGO(data, 'org.Hs.eg.db', ont = 'BP', pvalueCutoff = 0.05, pAdjustMethod = 'BY')
        BPdat = enrich_GO_BP.slots['result']
        
        robjects.r('k <- keys(org.Hs.eg.db, "ENTREZID")')
        robjects.r('df <- select(org.Hs.eg.db, keys=k, columns=c("ONTOLOGY","GO"),keytype="ENTREZID")')
        BP = robjects.r('df[which(df$ONTOLOGY == "BP"), ]')
        full_gl_GOBP = set(BP[0])
        full_gl_GOBP = set(map(int, full_gl_GOBP))
        
        df_GO_BP = pd.DataFrame(index=range(len(BPdat[0])))
        df_GO_BP['ID'] = BPdat[0]
        df_GO_BP['Description'] = BPdat[1]
        df_GO_BP['GeneRatio'] = BPdat[2]
        df_GO_BP['pvalue'] = BPdat[4]
        df_GO_BP['padjust'] = BPdat[5]
        df_GO_BP['Count'] = BPdat[8]
        
        empty_lijst_GO_BP = []
        for x in BPdat[7]: 
            lijst = x.split('/')
            lijst = list(map(int, lijst))
            empty_lijst_GO_BP.append(lijst)
            
        df_GO_BP['GeneID'] = empty_lijst_GO_BP
        df_GO_BP = df_GO_BP.sort_values("pvalue").head(500)
        
        dictionary = {}
        for index, row in df_GO_BP.iterrows():
            dictionary[row['ID']] = row['GeneID']
        
        return ProteinSet(dictionary, "GO_BP", full_gl_GOBP)
    
    def Reactome(self, clusterProfiler, ReactomePA, data): 
        """
        Function:
        ----------
        This function uses the R-package: ReactomePA to enrich a list of genes for Reactome pathways. You will get back a list of 
        Reactome pathways with the sets of genes that have something to do with that Reactome pathway. This list will be sorted 
        by the adjusted p-value. 
        
        Variable: 
        ----------
        ReactomePA = the R-package ReactomePA.
        data       = the gene list. 
        
        Returns: 
        ----------
        the call of the class: ProteinSet(). 
        """
        
        enrich_Reactome = ReactomePA.enrichPathway(gene=data, pvalueCutoff = 0.05, readable = True, pAdjustMethod = 'BY', organism = "human")
        ReactomeDat = enrich_Reactome.slots['result']
        
        reactome_db = rpackages.importr('reactome.db')
    
        reactome = robjects.r('as.list(reactomePATHID2EXTID)')
        list_ge = []
        list_ge = set(list_ge)
        for x in range(len(reactome)):
            new = [re.sub("\D", "", i) for i in list(reactome[x]) if i not in list_ge]
            list_ge = list_ge.union(new)
            
        full_gl_react = set(map(int, list_ge))
        
        df_Reactome = pd.DataFrame(index=range(len(ReactomeDat[0])))
        df_Reactome['ID'] = ReactomeDat[0]
        df_Reactome['Description'] = ReactomeDat[1]
        df_Reactome['GeneRatio'] = ReactomeDat[2]
        df_Reactome['pvalue'] = ReactomeDat[4]
        df_Reactome['padjust'] = ReactomeDat[5]
        df_Reactome['Count'] = ReactomeDat[8]
        
        entrez_list = []
        for x in ReactomeDat[7]:
            input_bitr = x.split('/')
            eg = clusterProfiler.bitr(input_bitr, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
            list_gene_string = list(eg[1])
            list_gene_int = list(map(int, list_gene_string))
            entrez_list.append(list_gene_int)
            
        df_Reactome['GeneID'] = entrez_list
        
        dictionary = {} 
        for index, row in df_Reactome.iterrows():
            dictionary[row['ID']] = row['GeneID']
        
        return ProteinSet(dictionary, "Reactome", full_gl_react)
    
    def make_dictio_ppi(self, get_ensp_filtered, get_entr_filtered_ens):
        """
        Function: 
        ----------
        This function maps the ppi dataset and the ageing related genes so that it can be enriched later on. 
        
        Variables: 
        ----------
        get_ensp_filtered = filtered ppi dataset turned into genes. 
        get_entr_filtered_ens = list of genes. 
        
        Returns: 
        ----------
        the call of the class: ProteinSet(). 
        """
        
        # Strip the 9606. from the identifiers!
        string = pd.read_csv('protein_links_v11.0_0.9.tsv', sep=' ')
        string['protein']  = string['protein'].map(lambda x: x.lstrip('9606.'))
        string['chemical'] = string['chemical'].map(lambda x: x.lstrip('9606.'))
        filt = string[string.protein.isin(get_ensp_filtered) | string.chemical.isin(get_ensp_filtered)]
        
        D = { p: [] for p in get_ensp_filtered }
        for i, r in filt.iterrows():
            if r.protein in get_ensp_filtered:
                D[r.protein]  = D[r.protein]  + [r.chemical]
            if r.chemical in get_ensp_filtered:
                D[r.chemical] = D[r.chemical] + [r.protein]
                
        D = { bm.protein_to_entrez(k) : [ bm.protein_to_entrez(q) for q in v] for (k,v) in D.items() }
        D = { k : [ e for e in v if e is not None] for (k,v) in D.items()}
        D = { k : v for (k,v) in D.items() if len(v) > 0 }
        
        return ProteinSet(D, "String", set(get_entr_filtered_ens))





