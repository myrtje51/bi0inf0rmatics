import sys
import pubchempy as pcp

def STITCH_inchikey():
    with open('protein_chemical_links_v5.0_2.0.tsv') as f2:
        next(f2)
        with_inchi = open("just_inchi.tsv","w+")
        with_inchi.write("InchiKey\n")
        for line2 in f2:
            fields = line2.strip().split()
            pubmedID_unfiltered = fields[0]
            pubmedID = pubmedID_unfiltered[4::]
            c = pcp.Compound.from_cid(pubmedID)
            inchikey = c.inchikey
            #print(inchikey)   
            with_inchi.write(inchikey + "\n")

        with_inchi.close() 
        
def main():
    STITCH_inchikey()

main() 
