'''
  Implementation of Neighbour joining to find a tree
  Using Progressive Pairwise dynamic programming to align sequences
  @author Simon Shrestha
  @version 8 December 2020
'''

import sys
from readfasta import readfasta
from tools import *
import time

# main menu for program
def printMenu( scoringList):
    print("\t\t NEIGHBOUR JOINING ")
    print(" A phylogenetic guide tree will be created using Neighbour Joining\n\n")
    print("-----------------------------------")
    print( "\t\t Scoring parameters: ")
    print( " Transitions score:     ", scoringList[0], "\t\tTransversions score:   ", scoringList[1])
    print( " Gap insertion penalty: ", scoringList[2], "\t\tGap Extension penalty: ", scoringList[3])
    print( " Match score:           ", scoringList[4])
    print("------------------------------------\n")
    print("\t\t Program Menu")
    print("\t\t______________")
    print(" 1. Start Program \t2. Edit Scoring Parameters")
    print(" 3. Exit")

def main():

    # Scoring parmaters for seq alignment score
    # transition, transversion, gap insertion, gap extension and match reward
    sl = [ -1, -2, -5, -3, 1]

    end = False
    while( not( end)):
        
        printMenu( sl)
        i = int( input(" Press the corresponding number\n "))

        # menu options
        
        if( i == 1):
            runCode = True
            inputFile = input(" Enter the input file of sequences in fasta format.\n :")
            try:
                inputlist = readfasta( inputFile)
            except FileNotFoundError:
                print( " invalid input filename or filepath")
                print( " Back to Main Menu \n\n")
                runCode = False

            if( runCode): 

                print( "\n")
                print( "\t Starting Program \n")
                print( " Input file:\t", inputFile)
                # to calculate time takes
                start = time.time()
                
                # prepare dna sequences list
                print(" Number of dna sequences:\t", len( inputlist))
                dnalist = []

                for i in inputlist:
                    dnalist.append( DNA( i[0], i[1]))

                # scoring list for distance measurement
                scoringList = [ -1, -2, -5, -3, 1]

                # make distance matrix
                m = matrixNJ( dnalist, scoringList)

                # feed the matrix to neighbour joining algorithm
                nj = neighbourJoining( m, dnalist)  

                end = time.time()
                print( " Time taken to calculate: ", end - start ," seconds")

            c = input("\n input 'y' to load program again any other characters to exit\n")
            if( c == 'y'):
                end = False
            else:
                end = True


        elif( i == 2):
            print( "scoring matrix")
            sl[0] = int( input( " Enter transition score\n "))
            sl[1] = int( input( " Enter transversion score\n "))
            sl[2] = int( input( " Enter gap insertion score\n "))
            sl[3] = int( input( " Enter gap extension score\n "))
            sl[4] = int( input( " Enter match reward score\n "))

        elif( i == 3):
            print( " Closing Program")
            end = True

        else:
            print( " Invalid input")
        
    
main()

    
