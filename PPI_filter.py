class PPIfilter:
    def get_PPIfilter(self):
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
    filtered = PPIfilter()
    print(filtered.get_PPIfilter()) 

main()
