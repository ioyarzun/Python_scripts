#!/usr/bin/python3

#exercise 1
class Protein(object):
    def __init__(self, identifier, sequence):
        self.identifier=str(identifier)
        self.sequence=str(sequence)

    def get_identifier(self):
        return self.identifier

    def get_sequence(self):
        return self.sequence

    def get_mw(self):
        aminoacids={'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16,'K': 146.19, 'M': 149.21,
        'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
        weight=0
        sequence=self.sequence.upper()
        for i in range(0,len(sequence)):
            weight=weight+aminoacids[sequence[i]]
        return format(weight, ".2f")

    def has_subsequence(self, subsequence):
        return subsequence.upper() in self.sequence

    def get_length(self):
        return len(self.sequence)
        

# test command
#my_prot=Protein("algo", "SDDSHDPISIGDSLDIGSHDHPSDYSILIDYD")
    #print(my_prot.get_identifier())

#exercise 2
def FASTA_iterator( fasta_filename ):
    fd=open(fasta_filename,"r")
    seq=""
    for line in fd:
        line=line.strip("\n")
        if line.startswith(">"):
            if seq:
                yield Protein(identifier,seq)
                seq=""
            identifier=line.strip(">")
        else:
            seq=seq+line
    yield Protein(identifier,seq)

# test command
#for protein in FASTA_iterator("example_fasta_file.fa.txt"):
    #print(protein.get_identifier())
