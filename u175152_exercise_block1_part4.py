#!/usr/bin/python3
#exercise 1
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


#exercise 2
def compare_fasta_file_identifiers( fasta_filenames_list ):
    frequency={}
    forspecific={}
    intersection=[]
    for fasta_file in fasta_filenames_list:
        # Calling the function from exercise number 1 in an object
        # that gives an identifier and its sequence each iteration
        forspecific[fasta_file]=[]
        actualfile=FASTA_iterator(fasta_file)
        for my_tuple in actualfile:
            identifier=my_tuple[0].upper()
            #Get the frequencies
            if identifier in frequency:
                frequency[identifier]+=1
            else:
                frequency[identifier]=1

            forspecific[fasta_file].append(identifier)

            if frequency[identifier]==len(fasta_filenames_list):
                intersection.append(identifier)

    union=set(frequency.keys())
    intersection=set(intersection)

    # List of identifiers that appears just in one file (freq=1)
    uniqueidentifiers=[]
    for identifier in frequency.keys():
        if frequency[identifier]==1:
            uniqueidentifiers.append(identifier)
    # Intersect the list with unique identifiers with the list of identifiers
    # on the different files
    specific={}
    for fasfile in forspecific:
        specific[fasfile]=set(forspecific[fasfile]).intersection(uniqueidentifiers)

    output_dictionary={"intersection":intersection, "union": union,
                        "frequency": frequency, "specific": specific}

    return output_dictionary

#test command
output=compare_fasta_file_identifiers(["example_fasta_file.fa.txt","fasta2.txt","fasta3.txt","fasta4.txt"])
print("this is output dictionary:      ", output)
