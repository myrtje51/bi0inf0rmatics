import sys
import pubchempy as pcp

def STITCH_inchikey():
    pubmedID_list = []
    with open('protein_chemical_links_v5.0_2.0.tsv') as f:
        next(f)
        for line in f:
            fields = line.strip().split()
            pubmedID_unfiltered = fields[0]
            pubmedID = pubmedID_unfiltered[4::]
            c = pcp.Compound.from_cid(pubmedID)
            print(c) 
            #pubmedID_list = []
            #pubmedID_list.append(pudmed_ID)
            #print(pubmedID_list) 

def main():
    STITCH_inchikey()

main() 
