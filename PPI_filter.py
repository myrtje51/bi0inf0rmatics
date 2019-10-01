class PPIfilter:
    def get_PPIfilter(self):
        """
        Function: 
        This function is supposed to give a file with all the protein-protein interactions that have a confidence-score higher than
        900 (0.9). This will be written to a different tsv file. 
        
        Variables:
        protein_links = the filte with all the protein-protein interactions (taken from STITCH) 
        ppis = the file where all the protein-protein interactions with a confidence-score of 900 (0.9) or higher will be written to. 
        ppi = the lines in the STRING file that is opened in the beginning. 
        fields_ppi = a way to get the same effect as awk has in unix. 
        combined_score_ppi = only the last column is selected from the line. Which happens to be the combined-score (also called the 
        confidence-score).
        """
        with open('9606.protein.links.v11.0.txt') as protein_links:
            next(protein_links)
            ppis = open("protein_links_v11.0_0.9.tsv","w+")
            ppis.write("chemical protein combined_score\n")
            for ppi in protein_links:
                fields_ppi = ppi.strip().split()
                combined_score_ppi = fields_ppi[2]
                if int(combined_score_ppi) > 899: 
                    ppis.write(ppi)
            ppis.close() 

def main():
    """
    Function: 
    to show how to call the class: PPIfilter. This class could also be called in another code, just don't forget to import the class. 
    """
    filtered = PPIfilter()
    print(filtered.get_PPIfilter()) 

main()
