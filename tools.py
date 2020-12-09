'''
  Methods to help find phylogenetic tree using
  neighbour joining  
  @author Simon Shrestha
  @version 7 December 2020
'''


# DNA validity checker
def isValid( seq):

    # valid characters
    validchar = ['A', 'C', 'G', 'T', '_', ' ']

    # change the sequence into Uppercase
    seq = seq.upper()

    # check every character in the sequence
    for c in seq:
        if (not (c in validchar)):
            return False

    # if no red flags
    return True

# global sequence alignment via dynamic programming
# firstSeq and secondSeq are the two sequences
# returns a list with two aligned sequences
def dpGlobal( firstSeq, secondSeq):

    # match reward, mismatch penalty and gap penalty
    match = 1
    mismatch = -1
    gap = -2

    # convert both sequences to all uppercase
    firstSeq = firstSeq.upper()
    secondSeq = secondSeq.upper()

    # m and n are the matrix dimensions
    m = len( firstSeq) + 1
    n = len( secondSeq) + 1

    # initialize mat as a matrix with value 0 in all positions
    mat = [ [0 for i in range(m)] for j in range(n)]

    # begin to fill the matrix
    for i in range(n):
        for j in range(m):
            # for position (0,0) set value as 0
            if j == 0 and i == 0:
                mat[i][j] = 0
            # set values of first row
            elif j == 0:
                mat[i][j] = mat[i - 1][j] + gap
            # set values of first column
            elif i == 0:
                mat[i][j] = mat[i][j - 1] + gap

            # for other positions calculate a value
            else:
                # set val as the mismatch penalty value
                val = mismatch
                # set value to match reward value
                # only if a match is found
                if( firstSeq[j - 1] == secondSeq[i - 1]):
                    val = match
                # set the maximum of the values from
                # added gap penalties from left or top
                # or added match or mismatch value
                mat[i][j] =  max( (mat[i][j-1] + gap),
                             (mat[i-1][j] + gap),
                             (mat[i-1][j-1] + val))

    # traceback

    # x and y are the variables to locate
    # positions in the matrix
    x = len( firstSeq)
    y = len( secondSeq)

    # initialize two strings to hold the
    # final sequence alignment
    firstSeqMod = ""
    secondSeqMod = ""

    # loops when both x and y are greater than 0
    while( x > 0 and y > 0):

        # these variables hold the calculated
        # values for traceback from
        # topleft, top and left positions
        val0 = mat[y][x] - mat[y - 1][x - 1]
        val1 = mat[y][x] - mat[y - 1][x]
        val2 = mat[y][x] - mat[y][x - 1]

        # when a path from top is found
        # add gap to the first sequence
        if( val1 == gap):
            y = y - 1
            firstSeqMod = "_" + firstSeqMod
            secondSeqMod = secondSeq[y] + secondSeqMod


        # when a path from left is found
        # add gap to the second seqence
        if( val2 == gap):
            x = x - 1
            firstSeqMod = firstSeq[x] + firstSeqMod
            secondSeqMod = "_" + secondSeqMod

        # when a path from topleft is found
        # add the matched or mismatched values
        if( val0 == match or val0 == mismatch):
            x = x - 1
            y = y - 1
            firstSeqMod = firstSeq[x] + firstSeqMod
            secondSeqMod = secondSeq[y] + secondSeqMod
           
    # until first column is reached
    # add gap to second sequence
    while( x > 0):
        x = x - 1
        firstSeqMod = firstSeq[x] + firstSeqMod
        secondSeqMod = "_" + secondSeqMod

    # until first row is reached
    # add gap to first sequence
    while( y > 0):
        y = y - 1
        firstSeqMod = "_" + firstSeqMod
        secondSeqMod = secondSeq[y] + secondSeqMod

    # prepare a return list
    # put the obtained sequences and return
    ret = []
    ret.append(firstSeqMod)
    ret.append(secondSeqMod)

    return ret

# create initial distance matrix to start NJ
# returns a distance matrix
def matrixNJ( dnaList, scoringList):

    # create a matrix of length of the dna list passed
    m = len( dnaList)
    mat = [ [0 for i in range(m)] for j in range(m)]

##    # scoring list for distance measurement
##    scoringList = [ -1, -2, -5, -3, 1]

    print( " Calculating pairwise distance\n")
    # calculate distance
    for i in range(m):
        for j in range(m):
            # calculate only top diagonal since the matrix will be symmetric
            if( i < j):
                mat[i][j] = DNA.dist( dnaList[i].seq, dnaList[j].seq, scoringList)
            if( j < i):
                mat[i][j] = mat[j][i]

    print( "Distance matrix obtained\n")
    for i in range(m):
        print ( mat[i])

    print( "--------------------------")
    return mat

# implementing neighbour joining
# takes distance matrix and list of dna
def neighbourJoining( dMatrix, dnaList):
    
    # create a copy of passed matrix
    localMat = dMatrix
    m = len( dnaList)
    r = 0

    print( " STARTING Neighbour Joining")
    print( "-----------------------------")

    print( " calculating divergence")

    # calculate divergence for each dna
    for i in range(m):
       for j in range(m):
           r = r + dMatrix[i][j]
           
       dnaList[i].setDivergence(r)
       r = 0

    # new distance is calculated for each pair
    for i in range(m):
        for j in range(m):
            if( i < j):
                localMat[i][j] = dMatrix[i][j] - (( dnaList[i].div + dnaList[j].div) / ( m - 2))
            if( j < i):
                localMat[i][j] = localMat[j][i]

    ret = []
    ret.append(reduceMatrix( localMat, dnaList))
    return ret

# join neighbouring sequences and reduce the matrix
# takes matrix and a data list as argument
def reduceMatrix( mat, dataList):

    # base case when there are only 2 nodes left in matrix
    if( len(mat) == 2):
        print(" complete")
        print(" --------------------")
        print(" Output")
        print(" --------------------")
        root = Node( dataList[0], dataList[1], 1, 1)
        root.name = "Root"
        toOutput( root)
        print(" \n\n")

        nlist =  preOrderTree( root)
        slist = ''.join( str(c) for c in nlist)

        print( slist)
        return root

    # s holds the smallest value initialized with a greater number
    # s = 1000000
    s = -100000
    m = len(mat)

    # find the location of the ##smallest largest distance value
    for i in range( m):
        for j in range( m):
            if( i < j):
                if (mat[i][j] > s):
                    s = mat[i][j]
                    a = i
                    b = j

    # save the distance for future use
    distAB = mat[a][b]

    # calculate internal distance from new Node
    distNodeA = (distAB) / 2 + (( dataList[a].div - dataList[b].div) / ( 2 * (m - 2)))
    distNodeB = (distAB) - distNodeA

    # remove the joined group              
    dnaA = dataList.pop(a)
    dnaB = dataList.pop(b - 1)

    # add the created node into data list
    dataList.append( Node( dnaA, dnaB, distNodeA, distNodeB))
    
##    print( "grouped:")
##    print( dnaA.name)
##    print( dnaB.name)
##    print( " -----")
    
    # save the distances from the joined node or dna
    distA = mat[a]
    distB = mat[b]

    # remove the values of merged dna or node
    for i in range(m):
        for j in range(m):
            if (j == a):
                mat[i].pop(j)

    for i in range(m):
        for j in range(m):
            if (j == ( b - 1)):
                mat[i].pop(j)
    mat.pop(a)
    mat.pop(b - 1)

    # new distance from new node
    nodeDist = []
    for i in range (len( distA)):
        newDist = ( distA[i] + distB[i] - distAB) / 2
        nodeDist.append( newDist)

    # append 0 at last since it is the distance to intself
    nodeDist.append( 0)

    # calculate the divegence of the node
    nodeDiv = 0
    for i in nodeDist:
        nodeDiv = nodeDiv + i

    dataList[ len(dataList) - 1].setDivergence( nodeDiv)

    # add node values in the distance matrix
    mat.append( nodeDist)

    for i in range (len( mat[0])):
        mat[i].append( mat[ len( mat) - 1][i]) 

    # recursive call for reducing matrix until only two nodes are left
    reduceMatrix( mat, dataList)

# method to traverse the tree in pre order
def preOrderTree( node):

    # if node is empty return
    if ( node == None):
        return ""

    # n list contains traversed nodes
    nlist = []
    nlist.append( node.name)
    nlist.append( "\n")

    pointer = "|___"

    hasChildB = False
    if( node.childB != None):
        hasChildB = True

    # traverse to children nodes
    visitNodes( node.childA, nlist, " ", pointer, hasChildB )
    visitNodes( node.childB, nlist, " ", pointer, False)

    return nlist

# traverse nodes in a tree
def visitNodes( node, nlist, padding, pointer, hasChildB):

    # if node is empty
    if( node != None):

        # stack necessary padding and name
        nlist.append( padding)
        nlist.append( pointer)
        nlist.append( node.name)
        nlist.append("\n")


        if (hasChildB):
            padding = padding + "|  "
        else:
            padding = padding + "  "
            
        pointer = "|___"

        hasChildB = False
        if( node.childB != None):
            hasChildB = True

        # recursive call to traverse further into tree
        visitNodes( node.childA, nlist, padding, pointer, hasChildB )
        visitNodes( node.childB, nlist, padding, pointer, False)

# output the obtained results
def toOutput( node):

    # information about grouped dna and the nodes connecting them
    if( node.type == "Node"):
        print(" NodeName = ", node.name)
        print("\tChildA = ", node.childA.name,"\n\tChildB = ", node.childB.name)
        toOutput( node.childA)
        toOutput( node.childB)    
    
# DNA class to create objects for holding dna informations
class DNA:

    # constructor
    def __init__( self, name, seq):
        self.type = "DNA"
        self.name = name
        if( isValid( seq)):
            self.seq = seq
        else:
            self.seq = "invalid sequence"
        self.childA = None
        self.childB = None

    # update the sequence
    def updateSeq( self, seq):
        if( isValid( seq)):
            self.seq = seq

    # set divergence
    def setDivergence( self, div):
        self.div = div

    # class method to find distance between two DNA sequences
    # with a passed scoring list to give score
    def dist( A, B, scoringList):

        # align the sequences using dynamic programming
        alignLst = dpGlobal( A, B)
        seqA = alignLst[0]
        seqB = alignLst[1]
        score = 0

        # scoring variables
        tranS = scoringList[0]
        tranV = scoringList[1]
        gapIns = scoringList[2]
        gapExt = scoringList[3]
        match = scoringList[4]

        # calculate the total score
        for i in range( len( seqA)):
            a = seqA[i]
            b = seqB[i]

            # match case
            if( a == b):
                score = score + match

            # gap found in first sequence
            elif( a == '_'):
                if( i == 0):
                    score = score + gapIns
                elif( seqA[ i - 1] == '_'):
                    score = score + gapExt
                else:
                    score = score + gapIns

            # gap found in second sequence
            elif( b == '_'):
                if( i == 0):
                    score = score + gapIns
                elif( seqB[ i - 1] == '_'):
                    score = score + gapExt
                else:
                    score = score + gapIns

            # transition penalties
            elif( (a == 'A' and b == 'G') or ( a == 'G' and b == 'A')):
                score = score + tranS

            elif( (a == 'C' and b == 'T') or ( a == 'T' and b == 'C')):
                score = score + tranS

            # transversion penalties
            else:
                score = score + tranV

        # return final score
        score = score
        return score

# Node class to hold node objects after joining two dna or nodes
class Node:

    # constructor
    def __init__( self, A, B, distToA, distToB):
        self.type = "Node"
        self.childA = A
        self.childB = B
        self.name = A.name[0] + B.name[0] + "Join"
        self.distToA = distToA
        self.distToB = distToB

    # set divergence for the node
    def setDivergence( self, div):
        self.div = div




        

        
        
        
        
        
    
    
