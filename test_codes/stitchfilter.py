class STITCH_filter: 
    def get_STITCH_filter(self):
        """
        Function:
        This function filters the human STITCH-data to only keeping the instances
        where the confidence score is 700 or higher. It then writes this data to
        a new file.
        
        Variables:
        file = the file with all the protein-chemical links that have a confidence-
        score of 700 or higher.
        f = the opened file with all the protein-chemical links in a human.
        fields = it is a way of imitating awk in python. This way you can select a
        column and do something with that column.
        combined_score = the last column a.k.a. the confidence score.
        """
        with open('9606.protein_chemical.links.v5.0.tsv') as f:
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
    """
    Calls the get_STITCH_filter() function. 
    """
    filtered_ST = STITCH_filter()
    print(filtered_ST.get_STITCH_filter()) 

main() 
