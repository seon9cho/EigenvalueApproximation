# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 15:28:51 2018

@author: Spencer
"""

import bestReductionClass as br
from matrixGenerator import matrixGenerator
import numpy as np
from time import time


#In every iteration, we generate a random matrix M and reduce along each row.
#For each row of M, we create a row of our training matrix X. Since our matrices
#have four rows, each iteration will produce four rows in X.

#Each row of the training matrix has four elements:
# 1) The flattened matrix M stored as a numpy array
# 2) The row we're reducing along (i.e. 0, 1, 2, or 3)
# 3) The original area of M
# 4) The new area of reduction of M

#The resulting matrix will have shape numIterations x 4 and is saved to the 
#local directory as "TrainingSet.npy" as a numpy array.

#Use the command 
# X = np.load("TrainingSet.npy") 
#to load the array back into the variable X.



#ADJUST THIS #########################################################################################
numIterations = 10000
######################################################################################################


variance1 = 4
variance2 = 1
matrixSize = 4

X = []
start = time()
print("Constructing Test Set...")
for i in range(numIterations):
    M = matrixGenerator(matrixSize, variance1, variance2)
    G1 = br.bestReductionClass(M)
    areas = G1.getAreas()
    for row in range(matrixSize):
        temp = np.array([M.flatten(), row, areas[0], areas[row + 1]])
        X.append(temp)
    if ((i + 1) % 500 == 0):
        print("Completed " + str(i + 1) + " iterations out of " + str(numIterations) +'.')
X = np.array(X)
np.save("TrainingSet.npy", X)
finish = time()
print("Finished")
print("It took " + str(finish - start) + " seconds to complete " + str(numIterations) + " iterations.")

