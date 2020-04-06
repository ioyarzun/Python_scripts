#!/usr/bin/python3

def FASTA_iterator( fasta_filename ):
    fd=open(fasta_filename,"r")
    seq=""
    for line in fd:
        line=line.strip("\n")
        if line.startswith(">"):
            if seq:
                yield (identifier,seq)
                seq=""
            identifier=line.strip(">")
        else:
            seq=seq+line
    yield (identifier,seq)

def get_max_sequence_length_from_FASTA_file (fasta_filename):
    generator=(len(sequence) for seq_id, sequence in FASTA_iterator(fasta_filename))
    return(max(generator))

print("Longest sequence length: ", get_max_sequence_length_from_FASTA_file("fasta2.txt"))

def get_min_sequence_length_from_FASTA_file ( fasta_filename ):
    generator=(len(sequence) for seq_id, sequence in FASTA_iterator(fasta_filename))
    return(min(generator))

print("Shortest sequence length: ", get_min_sequence_length_from_FASTA_file("fasta2.txt"))

def get_longest_sequences_from_FASTA_file( fasta_filename ):
    max_number=max(len(sequence) for seq_id, sequence in FASTA_iterator(fasta_filename))
    longestseq_tuples= list(sorted((seq_id.upper(),sequence) for seq_id, sequence in FASTA_iterator(fasta_filename) if len(sequence)==max_number))
    #longestseq_tuples.sort(key=lambda x: x[0]) using this we can remove sorted from previous line
    return(longestseq_tuples)

print("Longest sequence(s): ", get_longest_sequences_from_FASTA_file("fasta2.txt"))

def get_shortest_sequences_from_FASTA_file( fasta_filename ):
    min_number=min(len(sequence) for seq_id, sequence in FASTA_iterator(fasta_filename))
    shortestseq_tuples= list(sorted((seq_id.upper(),sequence) for seq_id, sequence in FASTA_iterator(fasta_filename) if len(sequence)==min_number))
    #shortestseq_tuples.sort(key=lambda x: x[0]) using this we can remove sorted from previous line
    return(shortestseq_tuples)

print("Shortest sequence(s): ",get_shortest_sequences_from_FASTA_file("fasta2.txt"))

def get_molecular_weights( fasta_filename ):
    aminoacids={'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16,'K': 146.19, 'M': 149.21,
    'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
    weight_dict={}
    for seq_id, sequence, seq_length in map(lambda x: (x[0].upper(), x[1], len(x[1])), ((seq_id, sequence) for seq_id, sequence in FASTA_iterator(fasta_filename))):
        weight=0
        weight_dict[seq_id]=0
        for i in range(0,seq_length):
            weight=weight+aminoacids[sequence[i]]
        weight_dict[seq_id]=format(weight, ".2f")
    return(weight_dict)

print("this is the weight dict ", get_molecular_weights("fasta2.txt"))

def get_sequence_with_max_molecular_weight( fasta_filename ):
    aminoacids={'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16,'K': 146.19, 'M': 149.21,
    'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
    weight=0
    for seq_id, sequence, seq_length in map(lambda x: (x[0].upper(), x[1], len(x[1])), ((seq_id, sequence) for seq_id, sequence in FASTA_iterator(fasta_filename))):
        temp_weight=0
        for i in range(0,seq_length):
            temp_weight= temp_weight+aminoacids[sequence[i]]
        if temp_weight>weight:
            weight=temp_weight
            max_weight_tuple=(seq_id,format(weight, ".2f"))
    return(max_weight_tuple)

print("this is the most weight tuple ", get_sequence_with_max_molecular_weight("fasta2.txt"))
