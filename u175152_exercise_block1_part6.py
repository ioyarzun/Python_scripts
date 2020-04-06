#!/usr/bin/python3

import sys
import math
import statistics

def calculate_pdb_chain_mean_minimum_distances(pdb_file_path=None):

    def PDB_iterator(pdb_file_path=None):
        if pdb_file_path==None:
            pdbfile=sys.stdin
        else:
            pdbfile=open(pdb_file_path,"r")
        for line in pdbfile:
            line=line.strip("\n")
            # Dividing fields in lines that start with ATOM
            if line.startswith("ATOM"):
                chain=line[21]
                residue=line[22:26]
                coordlist=[float(line[30:38].strip(" ")),float(line[38:46].strip(" ")),float(line[46:54].strip(" "))]
                yield(chain,residue,coordlist)
        pdbfile.close()

    my_dict={}
    for chain, resi, coordlist in  PDB_iterator(pdb_file_path):
        # Initialize chain dictionary inside my_dict
        if chain not in my_dict:
            my_dict[chain]={}
        if resi in my_dict[chain]:
            my_dict[chain][resi]=my_dict[chain][resi]+[coordlist]
        else: # Initialize resi list inside chain dictionary
            my_dict[chain][resi]=[coordlist]

    distance_dict={}
    for chain in my_dict:
        distance_dict[chain]=[]
        for resi in my_dict[chain]:
            for listofxyz in (my_dict[chain][ore] for ore in my_dict[chain].keys() if ore>resi):
                # List of xyz of any other residue in the chain
                # Important > for remove useless iteration between residues previously calculated
                distance=9999
                for x,y,z in my_dict[chain][resi]:
                    for ox,oy,oz in listofxyz:
                        temp_distance=math.sqrt((x-ox)**2+(y-oy)**2+(z-oz)**2)
                        if temp_distance<distance:
                            distance=temp_distance
                distance_dict[chain].append(distance)

    for chain in distance_dict:
        distance_dict[chain]=statistics.mean(distance_dict[chain])

    return(distance_dict)


# Calling the function defining the pdb_file in std input or with the
# pdb file path, if scipt is executed from the terminal
if __name__=='__main__':
    try:
        pdb_file=sys.argv[1]
    except:
        pdb_file=None
    diccionario=calculate_pdb_chain_mean_minimum_distances(pdb_file)
    for chain in diccionario:
        print(chain,":", format(diccionario[chain],".4f"))
