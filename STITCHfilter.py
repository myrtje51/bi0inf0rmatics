import sys

def STITCH_filter():
    with open('protein_chemical.links.v5.0.tsv') as f:
        next(f)
        file = open("protein_chemical_links_v5.0_2.1.tsv","w+")
        file.write("chemical\tprotein\tcombined_score\n")
        for line in f:
            fields = line.strip().split()
            combined_score = fields[2]
            if int(combined_score) > 699:
                file.write(line)
        file.close() 
def main():
    STITCH_filter()

main() 
