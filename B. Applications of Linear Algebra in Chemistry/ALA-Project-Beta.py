#Project Utility Kit
#Since this is a beta version, 
#some of the arguements in functions are explained in the other functions below them 
#Hence, bottom-up approach is better for traversal and understanding of functions 

import re
import math

#For taking transpose of matrix
def transpose(matrix):
    nRows = len(matrix)
    nCols = len(matrix[0])

    matrixT = list()

    for i in range(nCols):
        newRow = list()
        for j in range(nRows):
            newRow.append(matrix[j][i])
        matrixT.append(newRow)

    return matrixT

#For taking LCM of 2 numbers
def LCM(a,b):
   if a > b:
       greater = a
   else:
       greater = b

   while(True):
       if((greater % a == 0) and (greater % b == 0)):
           lcm = greater
           break
       greater += 1

   return lcm

#This function builds both the list of elements and adds the corresponding coefficients in the augmented matrix
def Add_Matrix(element, index, count, side):
    if index == len(eqnMatrix):
        eqnMatrix.append([])
        for i in elementsList:
            eqnMatrix[index].append(0)
    if element not in elementsList:
        elementsList.append(element)
        for j in range(len(eqnMatrix)):
            eqnMatrix[j].append(0)
    column = elementsList.index(element)
    eqnMatrix[index][column] += count*side

#This separates the finalCoeffns and elements to enable its addition in the respective matrices
def Element_Find(segment, index, atomsNum, side):
    elementQuantityNum = re.split('([A-Z][a-z]?)',segment) #storing numerical quantity of element
    i = 0
    while(i<len(elementQuantityNum)-1):
        i += 1
        if len(elementQuantityNum[i]) > 0:
            if elementQuantityNum[i+1].isdigit():
                count = int(elementQuantityNum[i+1])*atomsNum
                Add_Matrix(elementQuantityNum[i],index,count,side)
                i += 1
            else:
                Add_Matrix(elementQuantityNum[i],index,atomsNum,side)

#Following function splits the large and complex compounds 
#side indicates if it is a reactant or product referring to L.H.S. and R.H.S.
#For L.H.S. side = 1 ; for R.H.S. side = -1 
def DeComplex(compound,index,side):
    segments = re.split('(\([A-Za-z-0-9]*\)[0-9]*)',compound)
    for i in segments:
        #spliting elements inside a bracket
        if i.startswith("("):
            i = re.split('\)([0-9]*)',i)
            #storing the value of atoms for example Si(OH)4, OH has atomsNum value 4
            atomsNum = int(i[1])
            i = i[0][1:]
        else:
            atomsNum = 1
        Element_Find(i, index, atomsNum, side)
        


#----------------------------MAIN-------------------------
elementsList = [] #1D list that stores the individual distinct elements present in the unbalanced equation input by the user
eqnMatrix = [] #2D list that stores the augmented matrix for solving the equation


#Getting the unbalanced chemical equation as an input from the user, by taking two separate inputs for products and reactants.
print("Enter the LHS | Reactant Side \nNote: It is Case Sensitive | Leave No Spaces")
print("Example: H2+O2")
LHS = input()
print()
print("Enter the RHS | Products Side \nNote: It Is Case Sensitive | leave No Spaces")
print("Example: H2O")
RHS = input()

print()
print("Inputted Equation:")
print(LHS," --> ",RHS)
print()

# parting the reagents from each other and forming an array of reagents. For example: if the input is H2+O2 this gives an array LHS = ['H2','O2']
LHS = LHS.replace(' ','').split("+")
RHS = RHS.replace(' ','').split("+")


#Implementation

for i in range(len(LHS)):
    DeComplex(LHS[i],i,1)
for i in range(len(RHS)):
    DeComplex(RHS[i],i+len(LHS),-1)

#After this part we get the augmented matrix in eqnMatrix()
#Similarly List_element() would have the list of all the elements in the chemical reaction


#Here, we take the transpose of the finalCoeffns matrix inorder to get the correct augmented matrix, in which each row represents number of atoms of a element present in the equation
eqnMatrix = transpose(eqnMatrix)

from sympy import Matrix
eqnMatrix = Matrix(eqnMatrix)


# Solution finds the values of the coefficients in the balnced equation(Ax=0)
solution = eqnMatrix.nullspace()[0]
#print(solution)

#The process below checks if the solution is in fractional form, then to avoid fractional coefficients we multiply it with the appropriate lcm of the denominators
sol = []
for x in solution:
    sol.append(x)

num1 = sol[0].q
num2 = sol[1].q
#print(num1,num2)

find_lcm = LCM(num1,num2)
for i in range(2,len(sol)):
  find_lcm = LCM(find_lcm, sol[i].q)
multiplyFactor = find_lcm
solution = multiplyFactor*solution

#Generating the output
finalCoeffns = solution.tolist()
output = ""
for i in range(len(LHS)):
    output += str(finalCoeffns[i][0]) + LHS[i]
    if i<len(LHS) - 1:
        output += " + "
output += " --> "
for i in range(len(RHS)):
    output += str(finalCoeffns[i+len(LHS)][0]) + RHS[i]
    if i<len(RHS) - 1:
        output += " + "

print("Balanced Equation:")
print(output)