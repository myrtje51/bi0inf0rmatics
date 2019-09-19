import sys
import requests

def STITCH_inchikey():
    with open('protein_chemical_links_v5.0_2.1.tsv') as f2:
        next(f2)
        just_CIDs = open("just_CID.txt","w+")
        just_CIDs.write("cid=")
        count_line = 0
        for line2 in f2:
            count_line += 1
            fields = line2.strip().split()
            CID_unfiltered = fields[0]
            CID = CID_unfiltered[4::]
            if count_line < 466669 :
                just_CIDs.write(CID + ",")
            else:
                just_CIDs.write(CID)
        just_CIDs.close()
    data = open('just_CID.txt')
    response = requests.post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/InChIKey/CSV', data=data)
    CID_Inchi = open("resp_text.txt", "w+")
    CID_Inchi.write(response.text)
    CID_Inchi.close()

def main():
    STITCH_inchikey()

main()

