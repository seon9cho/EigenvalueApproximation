# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 11:05:23 2018

@author: Spencer
"""
import numpy as np
import sympy as sy
from sympy.printing.lambdarepr import LambdaPrinter
import time

class ImaginaryPrinter(LambdaPrinter):
    def _print_ImaginaryUnit(self, expr):
        return "1j"

class bestReduction():
    
    
    def __init__(self, matrix):
        self.x = sy.Symbol('x') 
        self.M = sy.Matrix(matrix)
        self.domain = self.getDomain(matrix)        
        
        
    def isospectralReduction(self, matrix, rowList):
        '''isospectral Reduction - isospectral reduction of matrix M over the
        over the index set rowList  
        
        @param rowList - a list specifying the rows of the matrix M to keep. E.g. to 
        keep rows 2, 3, and 4, we would input S = [1, 2, 3]. Recall that python
        indexing starts at zero.
        '''
        t1 = time.time()
        self.S = rowList #update S member
        S = list(rowList) 
        n = self.M.shape[0] #get dimension of M
        U = np.arange(0, n) #construct universal set
        SC = np.setxor1d(U, np.array(S)) # create complement array
                
        #get submatrices
        MSS = self.M[tuple(S), tuple(S),]
        MSSC = self.M[tuple(S), tuple(SC),]
        MSCSC = self.M[tuple(SC), tuple(SC),]
        MSCS = self.M[tuple(SC), tuple(S),]
        
        m = MSCSC.shape[0] #get height of MSCSC
        RSM = MSS - (1/(MSCSC[0]-self.x))*MSSC@MSCS
        t2 = time.time()
        print("Reduction: " + str(t2-t1))
        return RSM
    
    
    def lambdifyRegions(self, matrix):
        '''creates callable functions for our Gershgorwin Regions
        @returns m - A list of 2-tuples containing (f1, f2) where f1 and f2 are
        functions
        '''
        n = matrix.shape[0]
        m = []
        for i in range(n):
            f1 = sy.lambdify(self.x, abs(self.x - matrix.row(i)[i]))
            f2 = sy.lambdify(self.x, sum([abs(matrix.row(i)[j]) for j in range(n) if i != j]))
            m.append((f1,f2))            
        return m
    
    
    def estimateArea(self, matrix):
        ''' estimateArea - returns the number of points in our domain that are
        in our Gershgorwin Regions. The points in intersections are only 
        counted once.
        FIXME: Figure out if this can be vectorized somehow instead of using 
        for loops
        '''
        
        n = matrix.shape[0]
        t1 = time.time()
        m = self.lambdifyRegions(matrix)
        t2 = time.time()
        print("Regions: " + str(t2-t1))
        area = 0
        C = np.random.uniform(self.domain[0],self.domain[1],100) + 1j*np.random.uniform(self.domain[2],self.domain[3],100)
        for c in C:
            for i in range(n):
                if m[i][0](c) <= m[i][1](c): #if the point is in this region
                    area += 1
                    break #don't need to check the rest of the regions          
        return area
    
    
    def getDomain(self, matrix):
        '''
        getDomain - Automatically calculates the domain for Gershgorwin Regions
        and estimating area. 
        @paramter matrix
        @paramter density - an integer specifying the number of sample points
        per unit length
        '''
        n = matrix.shape[0]
        #get list of diagonals
        diags = matrix.diagonal()
        
        #get radii
        radii = []
        for i in range(n):
            radii.append(sum(abs(matrix[i][j]) for j in range(n) if i != j))
            
        
        maxDiagReal = max(diags.real)
        minDiagReal = min(diags.real)
        maxDiagImag = max(diags.imag)
        minDiagImag = min(diags.imag)
        maxRadius = max(radii)
        
        maxReal = maxDiagReal + maxRadius
        minReal = minDiagReal - maxRadius
        maxImag = maxDiagImag + maxRadius
        minImag = minDiagImag - maxRadius            

        return [minReal, maxReal, minImag, maxImag]
    
    
    def getPossibleSets(self, size=4):
        ''' returns a list of lists containing the correct indexes for 
        isospectral reductions
        '''
        U = np.arange(0, size)
        possibleRows = []
        for i in range(size):
            possibleRows.append(np.delete(U, i))
        return possibleRows
    
    
    def getAreas(self):
        t1 = time.time()
        n = self.M.shape[0]
        areas = []
        possibleRows = self.getPossibleSets(size=n)
        for S in possibleRows:
            RSM = self.isospectralReduction(self.M, S)
            areas.append(self.estimateArea(RSM))
        t2 = time.time()
        print("Calculating the areas: " + str(t2-t1))
        return areas            
    
    
    def syMatrixToNumpyArray(self, syMatrix):
        g = sy.lambdify((), syMatrix, modules="numpy", printer=ImaginaryPrinter)
        return g()
    
    
    
