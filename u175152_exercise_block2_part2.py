#!/usr/bin/python3
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
                raise ValueError("Impossible to create instance: " + letter + " not possible")

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
#

#Test commands
dna=DNASequence("asgdfs","atgatagatagaatagatatagatatagatagagatatagagatatagaacataagatatagagatagagaaagagctactagcaatgcatgcaatcatagactgacatgacatgcatgactgacaca")
twodna=dna+dna
rna=dna.transcribe()
#condna=dna.translate()
#conrna=rna.translate()

#print("________Ordered by molecular weight comprobation___________")
#list=[dna,rna,condna,conrna,dna]
#list.sort(reverse=True)
#for obj in list:
#    print(obj.get_mw(), "\t", obj.get_identifier())

# Adding objects in a diccionary as the key and random values
#print("__________Dictionary comprobation___________")
#q=ProteinSequence("ZZZ","atqmie")
#qd=ProteinSequence("ZZZ","atqmie")
#dicti={}
#dicti[q]='A'
#dicti[qd]='changed'
#dicti[conrna]='loquesea1'
#dicti[condna]='loquesea2'
#for object in dicti:
#    print(object.get_identifier(),dicti[object])

# Adding object to a set located in an object
#print("__________Set comprobation_________")
#prot=ProteinSequence("XXXX", "ahgtciemq")
#prot.add_object(q)
#prot.add_object(qd)
#prot.add_object(condna)
#prot.add_object(conrna)
#prot.add_object(prot)
#my_set=prot.objects_set
#len(my_set)
#for element in prot.objects_set:
#    print(element.get_identifier(), element.get_sequence())
