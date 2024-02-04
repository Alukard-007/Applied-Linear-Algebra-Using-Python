Project: Balancing Chemical Equations Using Concepts of Applied Linear Algebra in Python
Author : Buland Jayeshkumar Shah [Alukard.007]
Country: India
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

MAT248 Applied Linear Algebra | Project Report 
Applications of Linear Algebra in Chemistry | Balancing of Chemical Equations using Linear Algebraic Equations 


<Literature Review>
In the vast field of mathematics known as linear algebra, vectors, matrices, linear combinations, and vector spaces are all studied and implemented for efficiency. The topic of linear algebra has advanced, and it now finds use in many other contexts. In addition to the elimination method, matrices can be used to solve any system of linear equations. This idea can be used to balance chemical equations.

A chemical reaction is a process in which one or more compounds, known as reactants, are transformed into one or more distinct substances, known as products. Chemical equations use the symbols and formulas of the reactants and products to express a chemical reaction.

This project emphasizes the integration of knowledge from chemical equations and linear algebra to demonstrate one of its applications.


<Scope and Formulation of the Problem Statement>
One of the most fundamental and crucial ideas in chemistry is the concept of balancing chemical processes. Chemical reactions of various levels of complexity are present in every process around us. Chemical reactions are a part of everything, including breathing, digestion, and food preparation. These reactions involve compounds made up of different elements; some of these may be simple to manually balance, while others may be rather challenging. But in the end, it's just a straightforward linear equation that can be easily resolved using that method.

The law of conservation of matter, which is one of the fundamental laws, is what motivates the idea of balancing chemical equations and is presumed to hold true in all circumstances. The method of representing a chemical equation as a linear system is employed in the project to balance the equations.


<Methodology>
Based on research, the main focus will be on the role of linear algebra in balancing chemical equations with matrix notation. The project will hence demonstrate how to solve linear equations using matrices and matrix row operations. Given that it can handle any kind of chemical reaction involving a broad variety of reactants and products in this situation, the Gaussian Elimination approach comes out as an efficient algorithm for practical use in the project.

Following is the stepwise planned approach for our implementation, which is inspired on the basis of the research:

1. Taking the unbalanced chemical equation as input. 

>>> C2H5OH + O2 ---> CO2 + H20

2. Noting down the elements present in the compounds separately.

>>> C, H, O

3. Based on the elements, we will create column vectors for each compound, denoting the number of reacting elements in that compound. The number of rows in the column vector will be equal to the number of distinct elements that we noted in step 2. 

Similarly, the number of column vectors will depend on the number of compounds present on both sides of the reaction.

       C2H5OH   O2              CO2     H20
C:     | 2 |   | 0 |           | 1 |   | 0 |
H:     | 6 |   | 0 |           | 0 |   | 2 |
O:     | 1 |   | 2 |           | 2 |   | 1 |

Here, each column vector has the row 1, 2 and 3 denoting the numerical quantity of C, H and O respectively, for each compound. The first two representing that of the LHS side and the last two representing the RHS side.

4. Now, we assign coefficients for each compound on each side and assume the overall equation to be balanced. This is done in order to represent it in the parametric vector equation form.

(x1) C2H5OH + (x2) 02 ---> (x3) CO2 + (x4) H20

   | 2 |       | 0 |           | 1 |       | 0 |
x1 | 6 |  + x2 | 0 |  -->   x3 | 0 |  + x4 | 2 |
   | 1 |       | 2 |           | 2 |       | 1 |

5. We then equate the equation to zero and represent it in the form of an augmented matrix.

|  2  0  -1   0  0  |
|  6  0   0  -2  0  |
|  1  2  -2  -1  0  |

6. Using elimination techniques, we get to the Reduced Echelon Form of the matrix.

|    1      0      0   -(1/3)    0   |
|    0      1      0   -(1/1)    0   |
|    0      0      1   -(2/3)    0   |

7. From the above step, the values of coefficients are extracted and substituted in order to balance the chemical reaction.

8. The final balanced equation is sent as output.


<Motivation>
One of the most fundamental and crucial ideas in chemistry is the concept of balancing chemical reactions. Chemical reactions of various levels of complexity are present in every process around us. Chemical reactions are involved in respiration, digesting, and food preparation. These reactions involve compounds made up of different elements; some of these may be simple to manually balance, while others may be rather challenging. But in the end, it's just a straightforward linear equation that can be easily resolved using that method. Through the application of linear algebra knowledge, the initiative seeks to automate and streamline the process of balancing chemical reactions.


<Objectives>
The focus of this research will be on the role of linear algebra in balancing chemical equations with matrix notation. Additionally, it will show how to solve linear equations using matrices and matrix row operations. Also, it will be written in Python to save time and make it easier for everyone to use.

Given that it can handle any kind of chemical reaction involving a broad variety of reactants and products in this situation, the Gaussian Elimination approach is the best choice.
Through the use of computer technology, chemistry, and linear algebra, this project will be incredibly helpful in solving longer equations and will simplify the process overall.
For students who do average or even poorly, this offers a chance for progress. It may be able to eliminate what generally leads to failure and frustration for students who dislike chemistry. High achievers can also become incredibly quick and accurate, even in difficult situations. This study describes an augmented-matrix-based balancing technique.

Due to its unusual nature, it was best explained through examples in the process. Unquestionably, the matrix strategy has the benefit of being the most widely applicable technique for balancing chemical equations.


<References>
https://www.academia.edu/33623073/Balancing_of_Chemical_Equations_using_Matrix_Algebra 
https://www.wikihow.com/Balance-Chemical-Equations-Using-Linear-Algebra\ 
https://www.chemteam.info/Equations/Balance-Equation.html 
https://www.scirp.org/pdf/AM_2019071116274896.pdf

Output screenshots present in the directory are of the final project Python file.

Project was developed in between October - November 2022. 

~Sic Parvis Magna~
- Alukard.007
