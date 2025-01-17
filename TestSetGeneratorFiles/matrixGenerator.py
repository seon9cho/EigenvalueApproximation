# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 09:51:18 2018

@author: Spencer
"""

import numpy as np

def matrixGenerator(matrixSize = 4, variance1 = 4, variance2 = 1):
    '''
    Generates a random square matrix of size matrixSize. 
    
    The entries are generated as follows:
    For each diagonal entry Mii:
        mui = normal(0, variance1) (both real and complex)
        Mii = nomral(mui, variance1)
        
    For each row:
        mui' = normal(0, variance2)
        For each Mij, i != j:
            Mij = normal(mui', variance2)
    
    How much the Gershgorwin regions overlap depends on the relative sizes of 
    variance1 and variance 2. If variance1 >> variance2, we tend to get 
    nonoverlapping regions. If variance1 << variance2, we tend to get 
    overlapping regions. The default variances given allow us to generate a 
    good mix of matrices with overlapping and nonoverlapping regions.
    '''
    
    
    realDiagMeans = np.random.normal(0, variance1, 4)
    imagDiagMeans = np.random.normal(0, variance1, 4)
    realOffDiagMeans = np.random.normal(0, variance2, 4)
    imagOffDiagMeans = np.random.normal(0, variance2, 4)
    
    M = np.zeros((4, 4), dtype = complex)
    
    realDiags = []
    imagDiags = []
    for i in range(4):
        realDiags.append(np.random.normal(realDiagMeans[i], variance1))
        imagDiags.append(np.random.normal(imagDiagMeans[i], variance1))
    for i in range(matrixSize):
        M[i, :] = np.random.normal(realOffDiagMeans[i], variance2, 4) \
        +1j*np.random.normal(imagOffDiagMeans[i], variance2, 4)
    np.fill_diagonal(M, 0)
    M += np.diag(realDiags) + 1j*np.diag(imagDiags)
    return M
    