#!/usr/bin/python3
import sys
import os
import sequence_data as sd

def string_iterator(string,substrlength=1):
    i=0
    while i < len(string):
        if len(string[i:i+substrlength]) == substrlength:
            yield(string[i:i+substrlength])
        i+=substrlength

class Sequence(object):
    alphabet=""
    def __init__(self, identifier, sequence):
        self.__identifier=identifier.upper()
        self.__sequence=sequence.upper()
        self.objects_set=set() #initialize set for the object
        for letter in set(self.__sequence):
            if letter not in self.alphabet:
                raise(IncorrectSequenceLetter(letter, self.__class__.__name__))


    def get_identifier(self):
        return self.__identifier
    def get_sequence(self):
        return self.__sequence
    def get_mw(self):
        weight=0
        for i in range(0,len(self.__sequence)):
            weight=weight+self.weight_table[self.__sequence[i]]
        return float(format(weight, ".2f"))
    def has_subsequence(self, subsequence):
        return subsequence.upper() in self.__sequence

    def __len__(self):
        return len(self.__sequence)
    def __eq__(self, other_sequence):
        return  (self.__sequence == other_sequence.__sequence)
    def __ne__(self, other_sequence):
        return (self.__sequence != other_sequence.__sequence)
    def __add__(self,other):
        con_id=self.__identifier+"+"+other.__identifier
        con_seq=self.__sequence+other.__sequence
        return self.__class__(con_id,con_seq)
    def __getitem__(self,number):
        return self.__sequence[number]
    def __contains__(self,other):
        return other.__sequence in self.__sequence
    def __gt__(self,other):
        return self.get_mw() > other.get_mw()
    def __hash__(self):
        return hash((self.__identifier,self.__sequence))

    # Function for introduce objects in a set called objects previously
    # initialited. This set will be inside the "self" object.
    def add_object(self, obj):
        self.objects_set.add(obj)



class NucleotideSequence(Sequence):
    def string_iterator(string,substrlength=1):
        i=0
        while i < len(string):
            if len(string[i:i+substrlength]) == substrlength:
                yield(string[i:i+substrlength])
            i+=substrlength
    start_codons=[]
    stop_codons=[]
    def translate(self):
        translation_table={}
        trans_prot=""
        translation_activated=False
        for codon in string_iterator(self.get_sequence(), substrlength=3):
            if codon in self.start_codons:
                translation_activated=True
            if codon in self.stop_codons:
                translation_activated=False
            # Adding "" value for codons not defined in the translation table
            if codon not in self.translation_table.keys():
                translation_table[codon]=""
            if translation_activated:
                trans_prot=trans_prot+ self.translation_table[codon]
        return ProteinSequence(self.get_identifier(),trans_prot)
        # If self.__sequence or self.__identifier would be without __ we could call
        # sequence or identifier like objname.identifier o objname.sequence and not
        # with get_sequence() or get_identifier() functions
class DNASequence(NucleotideSequence):
    alphabet= 'GATC'
    start_codons=sd.dna_start_codons
    stop_codons=sd.dna_stop_codons
    translation_table=sd.dna_table
    weight_table=sd.dna_weights
    def transcribe(self):
        rna_seq=""
        dna_to_rna = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}
        for nucleotide in self.get_sequence():
            rna_seq=rna_seq+dna_to_rna[nucleotide]
        return RNASequence(self.get_identifier(), rna_seq)

class RNASequence(NucleotideSequence):
    alphabet = 'GAUC'
    start_codons=sd.rna_start_codons
    stop_codons=sd.rna_stop_codons
    translation_table=sd.rna_table
    weight_table=sd.rna_weights
    def reverse_transcribe(self):
        rna_to_dna = {'A': 'T', 'C': 'G', 'G': 'C', 'U': 'A'}
        dna_seq=""
        for nucleotide in self.get_sequence():
            dna_seq=dna_seq+rna_to_dna[nucleotide]
        return DNASequence(self.get_identifier(), dna_seq)

class ProteinSequence(Sequence):
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    weight_table=sd.protein_weights

class IncorrectSequenceLetter(ValueError):
    def __init__(self, letter, class_name):
        self.letter=letter
        self.class_name=class_name
    def __str__(self):
        return "The sequence item %s is not found in the alphabet of class %s" %(self.letter, self.class_name)

def FASTA_iterator( fasta_filename, class_name=ProteinSequence ):
    fd=open(fasta_filename,"r")
    seq=""
    for line in fd:
        line=line.strip("\n")
        if line.startswith(">"):
            if seq:
                try:
                    yield class_name(identifier,seq)
                    seq=""
                except IncorrectSequenceLetter as e:
                    sys.stderr.write(str(e)+"\n")
                    seq=""

            identifier=line.strip(">")
        else:
            seq=seq+line
    try:
        yield class_name(identifier,seq)
        seq=""
    except IncorrectSequenceLetter as e:
        sys.stderr.write(str(e)+"\n")
        seq=""

if __name__=='__main__':
    try:
        input=sys.argv[1]
        if input:
            try:
                output=sys.argv[2]
            except:
                output=None
    except:
        input=os.getcwd()
        output=None

    listofseqs=[]
    if os.path.isdir(input):
        os.chdir(input)
        i=sum(1 for file in os.listdir(input) if file.endswith((".fasta", ".fa")))
        sys.stderr.write(str(i)+" FASTA files found.\n")
        for file in os.listdir(input):
            if file.endswith((".fasta", ".fa")):
                for sequence in FASTA_iterator(file,class_name=DNASequence):
                    trans_seq=sequence.translate()
                    listofseqs.append(trans_seq)
                sys.stderr.write(input+file+" finished.\n")

    if os.path.isfile(input):
        i=0
        if input.endswith((".fasta", ".fa")):
            i+=1
        sys.stderr.write(str(i)+ " FASTA files found.\n")
        for sequence in FASTA_iterator(input,class_name=DNASequence):
            trans_seq=sequence.translate()
            listofseqs.append(trans_seq)

    sys.stderr.write(str(len(listofseqs)) + " sequences found.\n")
    sys.stderr.write("Sorting the sequences…\n")

    listofseqs.sort()

    sys.stderr.write("Sort process finished.\n")

    if output:
        output_file=open(output, "w")
        for seq in listofseqs:
            output_file.write("%s \t %s \t %s \n" %(seq.get_identifier(),len(seq),seq.get_mw()))
    else:
        for seq in listofseqs:
            print("%s \t %s \t %s \n" %(seq.get_identifier(),len(seq),seq.get_mw()))

    sys.stderr.write("Program finished correctly.\n")
    #ordenar por numero no string
