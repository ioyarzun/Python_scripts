
def count_sequences_by_residue_threshold(filename, residue, threshold=0.03):
    fd = open(filename,"r")
    protlength=0
    residuecount=0
    proteinnumber=0
    for line in fd:
        if line.startswith(">") and residuecount>0:
            residuefreq=residuecount/protlength
            protlength=0
            residuecount=0
            if residuefreq>threshold:
                proteinnumber+=1
        else:
            residuecount = residuecount + line.count(residue)
            protlength = protlength + len(line)
    residuefreq=residuecount/protlength #for the last protein
    if residuefreq>threshold:
        proteinnumber+=1
    return proteinnumber

#test command
print(count_sequences_by_residue_threshold("example_fasta_file.fa.txt","A", 0.03))


def print_sequence_tails(filename, output_filename, first_n=10, last_m=10):
#Get the proteins number
    fd = open(filename, "r")
    outfile=open(output_filename, "w")
    proteinnumber=0
    for line in fd:
        if line.startswith(">"):
            proteinnumber+=1
    outfile.write("#The file " + str(filename) + " contains " + str(proteinnumber) + " proteins. Here we show the code of the protein, the first " + str(first_n) + " aminoacids of each protein and the last " + str(last_m) + " aminoacids. \n\n")

#Get the Aa counts
    protseq=""
    print=""
    dd = open(filename, "r")
    for line in dd:
        line=line.strip("\n")
        if line.startswith(">") and len(protseq)>1:
            #Get and print the counts of the previous lines protein
            firstAa=protseq[0:first_n]
            lastAa=protseq[-last_m:]
            Aacounting=set(firstAa+lastAa) #Removing repeated Aa in the first and last Aas
            for Aa in Aacounting:
                Aacount=Aa+":"+str(protseq.count(Aa))
                print=print+ "," +Aacount

            outfile.write(protname+"\t"+firstAa+"\t"+lastAa+"\t"+print[1:]+"\n") #[1:] for delete the first ","

            protname=line.lstrip(">") #Get the name of the new protein
            protseq=""#reset previous protein seq
            print=""
        elif line.startswith(">"): #first protein name
            protname=line.lstrip(">")
        else:
            protseq=protseq+line

    #for the last protein
    firstAa=protseq[0:first_n]
    lastAa=protseq[-last_m:]
    Aacounting=set(firstAa+lastAa) #Removing repeated Aa in the first and last Aas
    for Aa in Aacounting:
        Aacount=Aa+":"+str(protseq.count(Aa))
        print=print+ "," +Aacount
    outfile.write(protname+"\t"+firstAa+"\t"+lastAa+"\t"+print[1:])
    outfile.close()


#test command
print_sequence_tails("example_fasta_file.fa.txt", "oyo.txt",3,6)
