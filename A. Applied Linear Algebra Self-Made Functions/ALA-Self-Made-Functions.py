"""
MAT248 : Applied Linear Algebra
Buland Shah [AU2140088]
"""

#__________________________
#Elementary Row Operations|______________________________________________________________________________

def interchange(matrix,row1Index,row2Index): #matrix is a 2D list with every element as a Row 
    matrix[row1Index], matrix[row2Index] = matrix[row2Index], matrix[row1Index]
    return matrix

def scaleBy(matrix,rowIndex,scalar):
    #for i in matrix[rowIndex]:
        #i = i*scalar
    for i in range(len(matrix[rowIndex])):
        matrix[rowIndex][i] = int(matrix[rowIndex][i] * scalar)
    return matrix

def addScaled(matrix,row1Index,row2Index,scalar):
    for i in range(len(matrix[row1Index])):
        matrix[row1Index][i] = matrix[row1Index][i] + (matrix[row2Index][i]*scalar)
    return matrix

#_________________________________________________________________________________________________________

#_____________
#Utility Kit |____________________________
#The Kit Contains Other Utility Functions|________________________________________________________________

def zeroCheck(matrix): #checks the zero rows and pushes them to the bottom of matrix
    nCols = len(matrix[0])
    nonzeroRows = list()
    leadingzeroRows = list() #to contain rows with leading enrty as zero
    for row in range(len(matrix)):
        RowElems = list(set(matrix[row]))
        if(len(RowElems)==1 and RowElems[0]==0):
            continue
        elif(matrix[row][0]==0):
            leadingzeroRows.append(matrix[row])
        else:
            #print(matrix[row])
            nonzeroRows.append(matrix[row])

    nonzeroRows.extend(leadingzeroRows)

    zeroRow = [0]*nCols
    numZeroRow = len(matrix) - len(nonzeroRows) #number of zero rows in matrix
    for i in range(numZeroRow):
        nonzeroRows.append(zeroRow)
    
    return nonzeroRows

def printMatrix(matrix):
    print()
    for row in range(len(matrix)):
        for col in range(len(matrix[0])):
            print(round(matrix[row][col]), end = "  ")
        print()

def zeroColCount(matrix): #checks for the presence of zero rows and returns list of zero column positions
    zeroCols = list()
    for j in range(len(matrix[0])):
        s=0
        for i in range(len(matrix)):
            s = s + abs(matrix[i][j])
        if(s==0):
            zeroCols.append(j)
    return zeroCols

def scaledRowEliminater(matrix): #checks for row that is equal to another scaled row and makes it zero row
    #for loop for taking rows into account 
    for i in range(len(matrix)):

        for j in range(i+1,len(matrix),1):

            count = 0 #counter for every row that counts the similar values between 2 rows (1st row--> i ; 2nd row--> j)
            for k in range(len(matrix[0])):
                if min(matrix[i][k], matrix[j][k])!= 0: #to avoid zero division error (x/0 is not defined)
                    scale = max(matrix[i][k], matrix[j][k]) / min(matrix[j][k], matrix[i][k])
                    if(scale - int(scale) == 0):
                        count = count + 1
            #print(count)
            if(count==len(matrix[0])):
                atleastOneScaled = True
                reqdRow = 0 #initialisation
                if(max(matrix[i][0], matrix[j][0])==matrix[i][0]):
                    reqdRow = i
                else:
                    reqdRow = j
                
                for z in range(len(matrix[0])):
                    matrix[reqdRow][z] = 0

    #printMatrix(matrix)
    matrix = zeroCheck(matrix)
    #printMatrix(matrix)
    return matrix

def diagonalOnes(matrix): #checks if its a matrix with diagonals 1 and others 0 (not necessarily a square matrix)
    diagOnes = 0
    Zeroes = 0
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if(i==j and abs(round(matrix[i][j]))==1):
                diagOnes = diagOnes + 1
            elif(i!=j and abs(round(matrix[i][j]))==0):
                Zeroes = Zeroes + 1

    totalElems = len(matrix) * len(matrix[0])

    return (diagOnes == len(matrix) and Zeroes == totalElems - diagOnes)

def smallMatrix(matrix, row, col): #to find minors for determinant calculation
    summatrix = [i[:col] + i[col+1:] for i in (matrix[:row] + matrix[row+1:])]
    return summatrix

def determinant(Matrix):
    #determinant only possible for square matrices
    if (len(Matrix) == len(Matrix[0])):  

        NumberOfRows = len(Matrix)
        for j in range(NumberOfRows):
            if NumberOfRows != len(Matrix[j]):
                print("Invalid matrix !")
                return 0       
        
        #recusrion ends at 2x2 matrix to get the minor:
        if(NumberOfRows == 2): 
            delta = (Matrix[0][0]*Matrix[1][1]) - (Matrix[0][1]*Matrix[1][0])
            return delta
        
        else:
            delta = 0
            for i in range(NumberOfRows):
                delta += ((-1)*i)*(Matrix[0][i])*(determinant(smallMatrix(Matrix,0,i)))

            return delta
    
    else:
        print("Invalid matrix !")
        return 0

#_________________________________________________________________________________________________________


#________________
#Major Functions|_________________________________________________________________________________________

#reduced=True returns the reduced echelon form
#solveEqn=True treats matrix as augmented matrix and neglects last col for consideration as pivot col
def toEchelon(matrix, showSteps=True, reduced=False, solveEqn=False, returnRank=False): 
    nRows = len(matrix)
    if(solveEqn):
        nCols = len(matrix[0])-1 #to neglect the last row
    else:
        nCols = len(matrix[0])

    if(showSteps):
        print('Original Matrix:')
        printMatrix(matrix)

    #pushing zero rows (if any) to the bottom
    matrix = zeroCheck(matrix)
    if(showSteps):
        print()
        print("Reaaranging Matrix:")
        printMatrix(matrix)

    #bringing the row with smallest leading entry to topmost
    minpos = 0
    c = 0 #minpos will hold the position of row with smallest leading entry
    for i in range(1,nRows,1):
        if(matrix[i][0]==0):
            continue
        elif(abs(matrix[i][0])<abs(matrix[minpos][0]) and abs(matrix[i][0])!=0):
            minpos = i
            c+=1
    if(c>0):
        matrix = interchange(matrix,minpos,0)

        if(showSteps):
            print()
            print("Row 1 <--> Row",minpos+1)
            printMatrix(matrix)

    #starting to find possible pivots
    #Basic Logic:
    #Take any non-zero col and search for non-zero value from row A to bottom row
    #Here, row A is the row present below the row having previous pivot or it is the first row
    #if we are looking for first pivot

    pivotsR = list() #stores row position of pivots
    pivotsC = list() #stores col position of pivots

    zColC = zeroColCount(matrix)
    for col in range(nCols):
        if(col not in zColC and col not in pivotsC):
            for row in range(nRows):
                if(row not in pivotsR and col not in pivotsC):
                    if(matrix[row][col]!=0):
                        #prevPivotR = pivotsR[-1]
                        if(len(pivotsR)!=0 and row!=pivotsR[-1]+1):
                            matrix = interchange(matrix,row,pivotsR[-1]+1)
                            if(showSteps):
                                print()
                                print("Row",row+1,"<--> Row",pivotsR[-1]+2)
                                printMatrix(matrix)
                            pivotsR.append(pivotsR[-1]+1)
                        else:
                            pivotsR.append(row)
                        pivotsC.append(col)

        #making values below current pivot as zeroes
        currPivotRow = pivotsR[-1]
        currPivotCol = pivotsC[-1]
        for h in range(currPivotRow+1,nRows,1):
            if(matrix[h][currPivotCol]!=0):
                #step i: making nonzero entry as zero by addition of scaled row
                #for simplicity we will make use of the row with the pivot
                reqdScalar = (-1) * (matrix[h][currPivotCol]) / (matrix[currPivotRow][currPivotCol])
                matrix = addScaled(matrix,h,currPivotRow,reqdScalar)
                if(showSteps):
                    print()
                    print("Row",h+1,"--> Row",h+1,"+ (",reqdScalar,")x Row",currPivotRow+1)
                    printMatrix(matrix)
            '''
            tmp = matrix    
            matrix = scaledRowEliminater(matrix)
            if(showSteps and tmp!=matrix):
                print("Deducting Similar/Scaled Rows:")
                printMatrix(matrix)
            '''


    #print(pivotsR)
    #print(pivotsC)
    #printMatrix(matrix)

    #printMatrix(matrix)

    #-----------------
    #The matrix we get here after all above operations is the one in echelon form
    #Now proceeding towards Reduced Echelon Form

    #Basic Logic
    #We take pivot columns and make every pivot 1 and 
    #then we we make all entries above pivot 1s as zeroes

    #step.i.
    #making all pivots as 1
    
    print()
    #print("REF")
    if(reduced):
        if(showSteps):
            print("Further Reducing to Reduced Echelon Form:\n")
        for k in range(len(pivotsC)):  
            if(matrix[pivotsR[k]][pivotsC[k]]!=1):
                scaleFactor = 1 / matrix[pivotsR[k]][pivotsC[k]]
                #print(matrix[pivotsR[k]][pivotsC[k]],scaleFactor)
                matrix = scaleBy(matrix,pivotsR[k],scaleFactor)
            if(showSteps):
                print()
                print("Row",pivotsR[k]+1,"= ( 1 /",matrix[pivotsR[k]][pivotsC[k]],") x Row",pivotsR[k]+1)
                printMatrix(matrix)

        for k in range(len(pivotsC)):
            
            #loop for moving up a pivot col from pivot to first row
            for i in range(pivotsR[k]-1,-1,-1):

                if(matrix[i][pivotsC[k]]!=0):
                    reqdScalar = (-1) * (matrix[i][pivotsC[k]]) // (matrix[pivotsR[k]][pivotsC[k]])
                    matrix = addScaled(matrix,i,pivotsR[k],reqdScalar)
                    if(showSteps):
                        print()
                        print("Row",i+1,"--> Row",i+1,"+ (",reqdScalar,")x Row",pivotsR[k]+1)
                        printMatrix(matrix)
                    
                    #promoting similar rows as zero rows
                    '''
                    tmp = matrix    
                    matrix = scaledRowEliminater(matrix)
                    if(showSteps and tmp!=matrix):
                        print("Deducting Similar/Scaled Rows:")
                        printMatrix(matrix)
                    '''

    if(returnRank):
        return len(pivotsR)
    else:
        return matrix
    
#------------------------------------------------------------------------------------------------------

def linearEqnSolver(matrix):
    import numpy as np
    #Dividing the Augmented Matrix in form of 
    #Ax = b
    Ax = list()
    b = list()

    for row in matrix:
        Ax.append(row[:-1])
        b.append(row[-1])
    #print(Ax)
    #print(b)

    if(diagonalOnes(Ax)):

        if(len(Ax)==len(Ax[0]) and np.linalg.det(Ax)!=0): #Matrix is Identity Matrix and det !=0
            #Unique Solution Exists
            cont = 1
            summ = 0 #to check if all unkowns are 0 or not
            print()
            for i in b:
                print("x",cont,"=",i)
                cont = cont + 1
                summ = summ + abs(i)
            print()
            if(summ == 0):
                print("This is a Trivial Solution.")
            else:
                print("This is a Non-Trivial Solution.")
        else:
            #Free Variable Exists
            print()
            for i in range(len(Ax)):
                print("x",i+1,"=",b[i])
            print()
            print("Other Unknown Variables Are Free")
            print("Hence, System of Linear Equations has more than one solution.")
    
    else:
        #Either No Soln (Inconsistent) or Infinitely Many Solutions
        inconsistent = False
        #infiniteSoln = False

        for i in range(len(Ax)):
            row = Ax[i]
            uniqueRowValues = list(set(row))
            if(len(uniqueRowValues)==1 and uniqueRowValues[0]==0):

                if(b[i]!=0):
                    inconsistent = True
                    break
        
        if(inconsistent):
            print("The System of Linear Equations is Inconsistent !")
        else:
            print("The System of Linear Equations has Infinitely More Solutions.")

#------------------------------------------------------------------------------------------------------

def norm(matrix):

    matrix = scaledRowEliminater(matrix)
    nRows = len(matrix)
    nCols = len(matrix[0])

    summ = 0
    for i in range(nRows):
        for j in range(nCols):
            summ = summ + ((matrix[i][j])**2)
    
    return round(summ**(0.5),2)

#------------------------------------------------------------------------------------------------------

def rank(matrix):
    return toEchelon(matrix, showSteps=False, reduced=True, returnRank=True)

#--------------------------------------------------------------------------------------------------------

def multiply(lhsMatrix, rhsMatrix):
    m1rows = len(lhsMatrix)
    m1cols = len(lhsMatrix[0])

    m2rows = len(rhsMatrix)
    m2cols = len(rhsMatrix[0])

    if(m1cols==len(rhsMatrix)):
        
            newMatrix = list()

            for k in range(m1rows):
                newMatrix.append([0]*len(rhsMatrix[0]))

            
            for i in range(m1rows):
                curr_row = lhsMatrix[i]

                for j in range(len(rhsMatrix[0])):
                    
                    summ = 0
                    for k in range(m1cols):
                        summ = summ + (curr_row[k]*rhsMatrix[k][j])
                    newMatrix[i][j] = summ

            return newMatrix

#---------------------------------------------------------------------------------------------------------

def inverse(Matrix):
    delta = determinant(Matrix)

    if(delta != 0): #condition for invertible matrix
        delta = abs(delta)
        nRows = len(Matrix)
        nCols = len(Matrix[0])
        for i in range(nRows):

            if nCols != len(Matrix[i]):
                print("Invalid Matrix !")
                return -1

        inverseMatrix = [[[0] for j in range(nCols)] for i in range(nRows)] #initialisation

        #calculating the cofactors of at ith,jth position and transposing with jth,ith position :-
        for i in range(nRows):

            for j in range(i,nCols):

                inverseMatrix[i][j] = ((-1)*(i+j))*(float(determinant(smallMatrix(Matrix, i, j)))*(1/delta))
                inverseMatrix[j][i] = ((-1)*(i+j))*(float(determinant(smallMatrix(Matrix, j, i)))*(1/delta))

                if(i != j):
                    inverseMatrix[i][j], inverseMatrix[j][i] = inverseMatrix[j][i], inverseMatrix[i][j]

        return inverseMatrix
    
    else:
        print("This matrix is not invertible :")
        return -1

#--------------------------------------------------------------------------------------------------------
#vList is a list containing vectors; 
#each vector in vList is in form:
#[x,y,z]-->for 3D or [0,0,x,y]->for 2D
def vectorMap(vList):
    from matplotlib import pyplot as plt
    import random

    dimension = len(vList[0])//2
    max = 0
    min = 0
    for i in range(len(vList)):
      for j in range(len(vList[0])):
          if(max<vList[i][j]):
              max=vList[i][j]
                
    for i in range(len(vList)):
      for j in range(len(vList[0])):
          if(min>vList[i][j]):
              min=vList[i][j]  

    if(dimension == 2):
       ax = plt.axes()
       for vector in vList:
           ax.arrow(vector[0], vector[1], vector[2], vector[3], head_width=0.5, head_length=0.5)
       plt.xlim(min-2,max+2)
       plt.ylim(min-2,max+2)
       plt.axvline(0, c='black', ls='--')
       plt.axhline(0, c='black', ls='--')

       plt.show()
    
    elif(dimension == 1):
        

        fig= plt.figure()
        ax= plt.axes(projection ='3d')
        ax.set_xlim([min,max])
        ax.set_ylim([min,max])
        ax.set_zlim([min,max])

        start=[0,0,0]
        colours = ['r','g','b','y']
        for vector in vList:
            a = random.randrange(0,3)
            ax.quiver(start[0],start[1],start[2],vector[0],vector[1],vector[2], color = colours[a])

        ax.view_init(10,10)
        plt.show()


#---------------------------------------------------------------------------------------------------------

# Sample Execution For Presentation 

print("Press 1 for Matrix Operations:")
print("Press 2 for Graphical Representation of Vectors:")
n = int(input())

if(n==1):
    print()
    print("---------------MATRIX OPERATIONS-----------------\n")
    #Matrix for Inverse and Row Reduction Operations
    a = [[0,1,0,2],[22,0,3,4],[0,0,0,0],[1,3,3,2]]
    b = [[5,4,3,2,1],[1,2,3,4,5],[1,3,5,7,9],[2,4,6,8,0],[2,3,5,6,1]]

    #Four Matrices for Multiplication
    m1 = [[1,2,1,2],[2,4,2,4],[1,1,1,1],[2,3,4,5]]
    m2 = [[2,2,2,2],[1,1,2,2],[2,5,6,7],[1,1,1,1]]

    m3 = [[1,2,1,2,1],[2,4,2,4,1],[1,1,1,1,1],[2,3,4,5,1],[2,3,4,2,1]]
    m4 = [[2,2,2,2,2],[1,1,2,1,2],[2,5,6,7,3],[1,1,1,1,1],[4,3,2,1,2]]

    print("Press 1 for Row Reductions and Corresponding Attributes:")
    print("Press 2 for Matrix Multiplication:")
    print("Press 3 for Solving a System of Linear Equations:")
    print("Press 4 for Matrix Inversion:")
    opt = int(input())

    if(opt==1):
        print("Demonstration of Echelon Form/s and its Attributes:\n")
        print("Case 1: 4x4 Matrix\n")
        toEchelon(a, reduced=True)

        print()
        print("Rank of Matrix:", rank(a))
        print("Norm of Matrix:",norm(a))

        print("-------------------------------------------------------------------------")
        print("Case 2: 5x5 Matrix\n")
        toEchelon(b, reduced=True)

        print()
        print("Rank of Matrix:", rank(b))
        print("Norm of Matrix:",norm(b))

    elif(opt==2):
        print("Matrix Multplication:\n")
        print("Case 1: 4x4 Matrices\n")
        print("Matrix 1:")
        printMatrix(m1)
        print()
        print("Matrix 2:")
        printMatrix(m2)
        print()
        print("Matrix 1 x Matrix 2:")
        printMatrix(multiply(m1,m2))
        print()
        print("--------------------------------------------------------------------------------------")

        print("Case 2: 5x5 Matrices\n")
        print("Matrix 1:")
        printMatrix(m3)
        print()
        print("Matrix 2:")
        printMatrix(m4)
        print()
        print("Matrix 1 x Matrix 2:")
        printMatrix(multiply(m3,m4))
        print()
        print("--------------------------------------------------------------------------------------")

    elif(opt==3):
        print("Solving a System of Linear Equations:\n")
        abc1 = [[2, -3, 1,-1],[1, -1, 2,-3],[3, 1, -1,9]]
        abc = [[8,7,0,9],[0,7,3,4],[2,5,2,5],[-8,2,3,1]]
        print("Case 1: 4x4 Matrix\n")
        #abc=[[1,1,1,8],[1,-1,2,6],[3,5,-7,14]]
        abc = toEchelon(abc, reduced=True, solveEqn=True)
        linearEqnSolver(abc)

        print()
        print("-----------------------------------------------------------------------------------------")
        print("Case 2: 5x5 Matrix\n")
        b = toEchelon(b, reduced=True, solveEqn=True)
        linearEqnSolver(b)

        print()
        print("-----------------------------------------------------------------------------------------")
        print("Case 3: 3x4  Matrix\n")
        abc1 = toEchelon(abc1, reduced=True, solveEqn=True)
        linearEqnSolver(abc1)
    
    elif(opt==4):
        #abc = [[2, -3, 1],[1, -1, 2],[3, 1, -1]]
        print("Matrix Inversion:")
        print()
        print("Case 1: 4x4 Matrix\n")
        print("Matrix A:")
        printMatrix(a)
        print()
        print("Matrix A Inverse:")
        inv = inverse(a)
        if(inv!=-1):
            for i in range(len(inv)):
                for j in range(len(inv[0])):
                    print(round(inv[i][j],3), end = '  ')
                print()
        print()

        print("Case 2: 5x5 Matrix\n")
        print("Matrix A:")
        printMatrix(b)
        print()
        print("Matrix A Inverse:")
        inv = inverse(b)
        if(inv!=-1):
            for i in range(len(inv)):
                for j in range(len(inv[0])):
                    print(round(inv[i][j],3), end = '  ')
                print()
        print()
    else:
        print("Invalid Input !")

elif(n==2):
    print()
    print("---------------VECTOR OPERATIONS-----------------\n")
    print("Press 1 for 2D Vector Graphical Represntation:")
    print("Press 2 for 3D Vector Graphical Represntation:")
    opt = int(input())

    if(opt==1):
        _2DvectorList = [[0,0,-10,15],[0,0,-10,2],[0,0,3,8]]
        print("2D Vector Representation:")
        print()
        count = 1
        for i in _2DvectorList:
            print("Vector:",count)
            print("Coordinates (x,y) :",i[2],",",i[3])
            count = count + 1

        print("Generating Graphical Plot...")
        vectorMap(_2DvectorList)
    
    elif(opt==2):
        _3DvectorList = [[2,3,1],[0,5,4],[7,-8,12]]
        print("3D Vector Representation:")
        print()
        count = 1
        for i in _3DvectorList:
            print("Vector:",count)
            print("Coordinates (x,y,z) :",i[0],",",i[1],",",i[2])
        
        print()
        print("Generating Graphical Plot...")
        vectorMap(_3DvectorList)

    else:
        print("Invalid Input !")
else:
    print("Invalid Input !")
    








#___________________
#<<<TESTING UNIT>>>|______________________________________________________________________________________

a = [[0,1,0],[22,0,3],[0,0,0],[1,3,3]]
b = [[0,1,0],[22,2,0],[0,0,0],[1,3,0]]
#print(toEchelon(a,reduced=True))
#print(zeroColCount(b))

c = [[0,1,0],[2,4,8],[11,22,33],[4,8,16]]
#scaledRowEliminater(c)

d = [[0,0,0,0,1,2],[1,1,1,2,1,0],[2,2,2,1,1,1]]
#d = [[0,5,4,4],[7,1,5,5],[2,3,8,8]]
#vectorMap(d)

e = [[1,2],[2,4]]
f = [[2],[2]]
#printMatrix(multiply(e,f))

#linearEqnSolver(a)

#_________________________________________________________________________________________________
#OUTPUT SCREEN

#1. EF/REF
'''
Demonstration of Echelon Form/s and its Attributes:

Case 1: 4x4 Matrix

Original Matrix:

0  1  0  2  
22  0  3  4  
0  0  0  0  
1  3  3  2  

Reaaranging Matrix:

22  0  3  4  
1  3  3  2  
0  1  0  2  
0  0  0  0  

Row 1 <--> Row 2

1  3  3  2  
22  0  3  4  
0  1  0  2  
0  0  0  0  

Row 2 --> Row 2 + ( -22.0 )x Row 1

1  3  3  2  
0  -66  -63  -40  
0  1  0  2  
0  0  0  0  

Row 3 --> Row 3 + ( 0.015151515151515152 )x Row 2

1  3  3  2  
0  -66  -63  -40  
0  0  -1  1  
0  0  0  0  

Further Reducing to Reduced Echelon Form:


Row 1 = ( 1 / 1 ) x Row 1

1  3  3  2  
0  -66  -63  -40  
0  0  -1  1  
0  0  0  0  

Row 2 = ( 1 / 1 ) x Row 2

1  3  3  2  
0  1  0  0  
0  0  -1  1  
0  0  0  0  

Row 3 = ( 1 / 1 ) x Row 3

1  3  3  2  
0  1  0  0  
0  0  1  -1  
0  0  0  0  

Row 1 --> Row 1 + ( -3 )x Row 2

1  0  3  2  
0  1  0  0  
0  0  1  -1  
0  0  0  0  

Row 1 --> Row 1 + ( -3 )x Row 3

1  0  0  5  
0  1  0  0  
0  0  1  -1  
0  0  0  0  


Rank of Matrix: 3
Norm of Matrix: 5.39
-------------------------------------------------------------------------
Case 2: 5x5 Matrix

Original Matrix:

5  4  3  2  1  
1  2  3  4  5  
1  3  5  7  9  
2  4  6  8  0  
2  3  5  6  1  

Reaaranging Matrix:

5  4  3  2  1  
1  2  3  4  5  
1  3  5  7  9  
2  4  6  8  0  
2  3  5  6  1  

Row 1 <--> Row 2

1  2  3  4  5  
5  4  3  2  1  
1  3  5  7  9  
2  4  6  8  0  
2  3  5  6  1  

Row 2 --> Row 2 + ( -5.0 )x Row 1

1  2  3  4  5  
0  -6  -12  -18  -24  
1  3  5  7  9  
2  4  6  8  0  
2  3  5  6  1  

Row 3 --> Row 3 + ( -1.0 )x Row 1

1  2  3  4  5  
0  -6  -12  -18  -24  
0  1  2  3  4  
2  4  6  8  0  
2  3  5  6  1  

Row 4 --> Row 4 + ( -2.0 )x Row 1

1  2  3  4  5  
0  -6  -12  -18  -24  
0  1  2  3  4  
0  0  0  0  -10  
2  3  5  6  1  

Row 5 --> Row 5 + ( -2.0 )x Row 1

1  2  3  4  5  
0  -6  -12  -18  -24  
0  1  2  3  4  
0  0  0  0  -10  
0  -1  -1  -2  -9  

Row 3 --> Row 3 + ( 0.16666666666666666 )x Row 2

1  2  3  4  5  
0  -6  -12  -18  -24  
0  0  0  0  0  
0  0  0  0  -10  
0  -1  -1  -2  -9  

Row 5 --> Row 5 + ( -0.16666666666666666 )x Row 2

1  2  3  4  5  
0  -6  -12  -18  -24  
0  0  0  0  0  
0  0  0  0  -10  
0  0  1  1  -5  

Row 5 <--> Row 3

1  2  3  4  5  
0  -6  -12  -18  -24  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  

Further Reducing to Reduced Echelon Form:


Row 1 = ( 1 / 1 ) x Row 1

1  2  3  4  5  
0  -6  -12  -18  -24  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  

Row 2 = ( 1 / 1 ) x Row 2

1  2  3  4  5  
0  1  2  3  4  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  

Row 3 = ( 1 / 1.0 ) x Row 3

1  2  3  4  5  
0  1  2  3  4  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  

Row 4 = ( 1 / 1 ) x Row 4

1  2  3  4  5  
0  1  2  3  4  
0  0  1  1  -5  
0  0  0  0  1  
0  0  0  0  0  

Row 1 --> Row 1 + ( -2 )x Row 2

1  0  -1  -2  -3  
0  1  2  3  4  
0  0  1  1  -5  
0  0  0  0  1  
0  0  0  0  0  

Row 2 --> Row 2 + ( -2.0 )x Row 3

1  0  -1  -2  -3  
0  1  0  1  14  
0  0  1  1  -5  
0  0  0  0  1  
0  0  0  0  0  

Row 1 --> Row 1 + ( 1.0 )x Row 3

1  0  0  -1  -8  
0  1  0  1  14  
0  0  1  1  -5  
0  0  0  0  1  
0  0  0  0  0  

Row 3 --> Row 3 + ( 5.0 )x Row 4

1  0  0  -1  -8  
0  1  0  1  14  
0  0  1  1  0  
0  0  0  0  1  
0  0  0  0  0  

Row 2 --> Row 2 + ( -14.0 )x Row 4

1  0  0  -1  -8  
0  1  0  1  0  
0  0  1  1  0  
0  0  0  0  1  
0  0  0  0  0  

Row 1 --> Row 1 + ( 8.0 )x Row 4

1  0  0  -1  0  
0  1  0  1  0  
0  0  1  1  0  
0  0  0  0  1  
0  0  0  0  0  


Rank of Matrix: 4
Norm of Matrix: 2.65
'''

#-------------------------------------------------------------------------------------------------------
#2. Matrix Multiplication

'''
Matrix Multplication:

Case 1: 4x4 Matrices

Matrix 1:

1  2  1  2  
2  4  2  4  
1  1  1  1  
2  3  4  5  

Matrix 2:

2  2  2  2  
1  1  2  2  
2  5  6  7  
1  1  1  1  

Matrix 1 x Matrix 2:

8  11  14  15  
16  22  28  30  
6  9  11  12  
20  32  39  43  

--------------------------------------------------------------------------------------
Case 2: 5x5 Matrices

Matrix 1:

1  2  1  2  1  
2  4  2  4  1  
1  1  1  1  1  
2  3  4  5  1  
2  3  4  2  1  

Matrix 2:

2  2  2  2  2  
1  1  2  1  2  
2  5  6  7  3  
1  1  1  1  1  
4  3  2  1  2  

Matrix 1 x Matrix 2:

12  14  16  14  13  
20  25  30  27  24  
10  12  13  12  10  
24  35  41  41  29  
21  32  38  38  26  

--------------------------------------------------------------------------------------
'''

#---------------------------------------------------------------------------------------------
#3. Solving Linear Equations

'''
Solving a System of Linear Equations:

Case 1: 4x4 Matrix

Original Matrix:

8  7  0  9  
0  7  3  4  
2  5  2  5  
-8  2  3  1  

Reaaranging Matrix:

8  7  0  9  
2  5  2  5  
-8  2  3  1  
0  7  3  4  

Row 1 <--> Row 2

2  5  2  5  
8  7  0  9  
-8  2  3  1  
0  7  3  4  

Row 2 --> Row 2 + ( -4.0 )x Row 1

2  5  2  5  
0  -13  -8  -11  
-8  2  3  1  
0  7  3  4  

Row 3 --> Row 3 + ( 4.0 )x Row 1

2  5  2  5  
0  -13  -8  -11  
0  22  11  21  
0  7  3  4  

Row 3 --> Row 3 + ( 1.6923076923076923 )x Row 2

2  5  2  5  
0  -13  -8  -11  
0  0  -3  2  
0  7  3  4  

Row 4 --> Row 4 + ( 0.5384615384615384 )x Row 2

2  5  2  5  
0  -13  -8  -11  
0  0  -3  2  
0  0  -1  -2  

Row 4 --> Row 4 + ( -0.5151515151515151 )x Row 3

2  5  2  5  
0  -13  -8  -11  
0  0  -3  2  
0  0  0  -3  

Further Reducing to Reduced Echelon Form:


Row 1 = ( 1 / 1 ) x Row 1

1  2  1  2  
0  -13  -8  -11  
0  0  -3  2  
0  0  0  -3  

Row 2 = ( 1 / 1 ) x Row 2

1  2  1  2  
0  1  0  0  
0  0  -3  2  
0  0  0  -3  

Row 3 = ( 1 / 1 ) x Row 3

1  2  1  2  
0  1  0  0  
0  0  1  0  
0  0  0  -3  

Row 1 --> Row 1 + ( -2 )x Row 2

1  0  1  2  
0  1  0  0  
0  0  1  0  
0  0  0  -3  

Row 1 --> Row 1 + ( -1 )x Row 3

1  0  0  2  
0  1  0  0  
0  0  1  0  
0  0  0  -3  
The System of Linear Equations is Inconsistent !

-----------------------------------------------------------------------------------------
Case 2: 5x5 Matrix

Original Matrix:

5  4  3  2  1  
1  2  3  4  5  
1  3  5  7  9  
2  4  6  8  0  
2  3  5  6  1  

Reaaranging Matrix:

5  4  3  2  1  
1  2  3  4  5  
1  3  5  7  9  
2  4  6  8  0  
2  3  5  6  1  

Row 1 <--> Row 2

1  2  3  4  5  
5  4  3  2  1  
1  3  5  7  9  
2  4  6  8  0  
2  3  5  6  1  

Row 2 --> Row 2 + ( -5.0 )x Row 1

1  2  3  4  5  
0  -6  -12  -18  -24  
1  3  5  7  9  
2  4  6  8  0  
2  3  5  6  1  

Row 3 --> Row 3 + ( -1.0 )x Row 1

1  2  3  4  5  
0  -6  -12  -18  -24  
0  1  2  3  4  
2  4  6  8  0  
2  3  5  6  1  

Row 4 --> Row 4 + ( -2.0 )x Row 1

1  2  3  4  5  
0  -6  -12  -18  -24  
0  1  2  3  4  
0  0  0  0  -10  
2  3  5  6  1  

Row 5 --> Row 5 + ( -2.0 )x Row 1

1  2  3  4  5  
0  -6  -12  -18  -24  
0  1  2  3  4  
0  0  0  0  -10  
0  -1  -1  -2  -9  

Row 3 --> Row 3 + ( 0.16666666666666666 )x Row 2

1  2  3  4  5  
0  -6  -12  -18  -24  
0  0  0  0  0  
0  0  0  0  -10  
0  -1  -1  -2  -9  

Row 5 --> Row 5 + ( -0.16666666666666666 )x Row 2

1  2  3  4  5  
0  -6  -12  -18  -24  
0  0  0  0  0  
0  0  0  0  -10  
0  0  1  1  -5  

Row 5 <--> Row 3

1  2  3  4  5  
0  -6  -12  -18  -24  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  

Further Reducing to Reduced Echelon Form:


Row 1 = ( 1 / 1 ) x Row 1

1  2  3  4  5  
0  -6  -12  -18  -24  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  

Row 2 = ( 1 / 1 ) x Row 2

1  2  3  4  5  
0  1  2  3  4  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  

Row 3 = ( 1 / 1.0 ) x Row 3

1  2  3  4  5  
0  1  2  3  4  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  

Row 1 --> Row 1 + ( -2 )x Row 2

1  0  -1  -2  -3  
0  1  2  3  4  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  

Row 2 --> Row 2 + ( -2.0 )x Row 3

1  0  -1  -2  -3  
0  1  0  1  14  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  

Row 1 --> Row 1 + ( 1.0 )x Row 3

1  0  0  -1  -8  
0  1  0  1  14  
0  0  1  1  -5  
0  0  0  0  -10  
0  0  0  0  0  
The System of Linear Equations is Inconsistent !

-----------------------------------------------------------------------------------------
Case 3: 3x4  Matrix

Original Matrix:

2  -3  1  -1  
1  -1  2  -3  
3  1  -1  9  

Reaaranging Matrix:

2  -3  1  -1  
1  -1  2  -3  
3  1  -1  9  

Row 1 <--> Row 2

1  -1  2  -3  
2  -3  1  -1  
3  1  -1  9  

Row 2 --> Row 2 + ( -2.0 )x Row 1

1  -1  2  -3  
0  -1  -3  5  
3  1  -1  9  

Row 3 --> Row 3 + ( -3.0 )x Row 1

1  -1  2  -3  
0  -1  -3  5  
0  4  -7  18  

Row 3 --> Row 3 + ( 4.0 )x Row 2

1  -1  2  -3  
0  -1  -3  5  
0  0  -19  38  

Further Reducing to Reduced Echelon Form:


Row 1 = ( 1 / 1 ) x Row 1

1  -1  2  -3  
0  -1  -3  5  
0  0  -19  38  

Row 2 = ( 1 / 1 ) x Row 2

1  -1  2  -3  
0  1  3  -5  
0  0  -19  38  

Row 3 = ( 1 / 1 ) x Row 3

1  -1  2  -3  
0  1  3  -5  
0  0  1  -2  

Row 1 --> Row 1 + ( 1 )x Row 2

1  0  5  -8  
0  1  3  -5  
0  0  1  -2  

Row 2 --> Row 2 + ( -3 )x Row 3

1  0  5  -8  
0  1  0  1  
0  0  1  -2  

Row 1 --> Row 1 + ( -5 )x Row 3

1  0  0  2  
0  1  0  1  
0  0  1  -2  

x 1 = 2
x 2 = 1
x 3 = -2

This is a Non-Trivial Solution.
'''

#----------------------------------------------------------------------------------------------------
#4. Inverse of a Matrix 

'''
Matrix Inversion:

Case 1: 4x4 Matrix

Matrix A:

0  1  0  2  
22  0  3  4  
0  0  0  0  
1  3  3  2  

Matrix A Inverse:
This matrix is not invertible :

Case 2: 5x5 Matrix

Matrix A:

5  4  3  2  1  
1  2  3  4  5  
1  3  5  7  9  
2  4  6  8  0  
2  3  5  6  1  

Matrix A Inverse:
0.0  -0.024  -0.027  0.126  0.279  
0.067  0.051  0.043  0.253  0.28  
0.116  0.086  0.064  0.294  0.305  
0.113  0.088  0.061  0.318  0.328  
0.157  0.179  0.128  0.282  0.237  

'''

#--------------------------------------------------------------------------------------------------
