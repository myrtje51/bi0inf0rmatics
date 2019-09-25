import sys
import requests

def STITCH_inchikey():
    """
    Function:
    This function returns a list of CIDs from STITCH and their responsive
    InChiKeys. It first puts all the CIDs in a format that can be read by
    PubChem and then uses requests to imitate the "curl" command from Bash
    and get an output.
    
    Variables:
    just_CIDs = a file with just a list full of CIDs seperated by a comma.
    count_line = counts which line the for-loop is reading at this moment.
    f2 = the file with the filtered protein-chemical links.
    fields = this imitates the effect awk has in bash. This way you can select
    a column and do something with this "value"
    CID_unfiltered = the CID in the f2 file looks like this: CIDm000001. PubChem
    asks for an input that looks like this: 000001 (so without the CIDm). This
    needs to be filtered then.
    CID = the CID that is without the CIDm. This can be used as input for
    PubChem.
    data = the data that is going to be used in the PubChem curl.
    response = the response of the curl from PubChem in csv format.
    CID_inchi = a file with the output from response (which has the CID and the
    InChiKey in one table). 
    """
    with open('protein_chemical_links_v5.0_2.1.tsv') as f3:
        next(f3)
        size=sum(1 for _ in f3)
        size = size + 1 
    with open('protein_chemical_links_v5.0_2.1.tsv') as f2: 

        just_CIDs = open("just_CID.txt","w+")
        just_CIDs.write("cid=")
        count_line = 0
        for line2 in f2:

            count_line += 1 
            fields = line2.strip().split()
            CID_unfiltered = fields[0]
            CID = CID_unfiltered[4::]
            if count_line < size:

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
    """
    Function:
    Calls the STITCH_inchikey() function. 
    """
    STITCH_inchikey()

main()

