# NeighbourJoining

Program to find a phylogenetic tree using Neighbour Joining algorithm using progressive pairwise dynamic
programming to align sequences.

to start program:
  run NJ.py

  an interface will be shown aith various options.

input files:
  takes in DNA sequences in fasta format

scoring parameters:
  the default scoring parameters for
  transitions: -1 , transversions: -2,
  gap insertions: -5 gap extensions: -3,
  match reward: 1

  These prameters can be manually set through the interface.

required files:
  tools.py  contains the functions and classes to claculate and find phylogenetic tree.
  readfasta.py  contains the function to read file in fasta format.
  
output:
  prints the provided sequences and the calculated tree.
