#!/usr/bin/python3

def calculate_aminoacid_frequencies(fasta_filename,subsequences_filename,number_of_repetitions,output_filename):

#Get the sequences in an array and the proteins number from fasta file
    fd = open(fasta_filename,"r")
    proteinnumber=0
    arrayofsequences=[]
    for line in fd:
        line=line.strip("\n")
        if line.startswith(">"):
            proteinnumber+=1
            arrayofsequences.append("")

        else:
            arrayofsequences[proteinnumber-1]=arrayofsequences[proteinnumber-1]+line

#Get the number of sequences repeated in each protein in a dictionary
    sf = open(subsequences_filename, "r")
    my_dict= dict()
    subseqnumber=0
    for subsequence in sf:
        subsequence=subsequence.strip("\n")
        my_dict[str(subsequence)]=0
        subseqnumber+=1
        for protsequence in arrayofsequences:
            repetitions=protsequence.count(subsequence)
            if repetitions >= number_of_repetitions:
                my_dict[subsequence]+=1

#Print the dictionary reverse ordered by the values of the keys
    my_dictord = dict(sorted(my_dict.items(), reverse=True, key=lambda x: x[1]))
    outfile=open(output_filename, "w")
    outfile.write("#Number of proteins:\t\t{1:10d}{0}#Number of subsequences:\t{2:10d}{0}#Subsequence proportions:{0}".format("\n",proteinnumber,subseqnumber))
    for key in my_dictord:
        outfile.write("{1}\t{2:>10d}\t{3:.4f}{0}".format("\n",key,my_dictord[key],my_dictord[key]/proteinnumber))
    outfile.close()

#test command
calculate_aminoacid_frequencies("example_fasta_file.fa.txt","sequence_fragments.txt",2, "outputfile.txt")
