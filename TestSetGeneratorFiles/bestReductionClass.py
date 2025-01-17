# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 11:05:23 2018

@author: Spencer
"""
import numpy as np
import sympy as sy
import math
from sympy.printing.lambdarepr import LambdaPrinter
#import time

class ImaginaryPrinter(LambdaPrinter):
    def _print_ImaginaryUnit(self, expr):
        return "1j"

class bestReductionClass():
    
    
    def __init__(self, matrix):
        self.x = sy.Symbol('x') 
        self.M = sy.Matrix(matrix)
        self.numSamples = 100
        self.realRange = -1
        self.imagRange = -1
        self.domain = self.getDomain(matrix)           
        
        
    def isospectralReduction(self, matrix, rowList):
        '''isospectral Reduction - isospectral reduction of matrix M over the
        over the index set rowList  
        
        @param rowList - a list specifying the rows of the matrix M to keep. E.g. to 
        keep rows 2, 3, and 4, we would input S = [1, 2, 3]. Recall that python
        indexing starts at zero.
        '''
        #start = time.time()
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
        RSM = MSS - MSSC*(MSCSC - self.x*sy.eye(m)).inv()*MSCS
        #finish = time.time()
        #print("time for isospectral reduction: " + str(finish - start))
        return RSM
    
    
    def lambdifyRegions(self, matrix):
        '''creates callable functions for our Gershgorwin Regions
        @returns m - A list of 2-tuples containing (f1, f2) where f1 and f2 are
        functions
        '''
        #start = time.time()
        n = matrix.shape[0]
        m = []
        for i in range(n):
            f1 = sy.lambdify(self.x, abs(self.x - matrix.row(i)[i]))
            f2 = sy.lambdify(self.x, sum([abs(matrix.row(i)[j]) for j in range(n) if i != j]))
            m.append((f1,f2))   
        #end = time.time()
        #print("time for lamdifyRegions: " + str(end - start))
        return m
    
    
    def estimateArea(self, matrix):
        ''' estimateArea - returns the geometric area of the intersection of 
        the Gershgorwin Regions of the given matrix. The points in 
        intersections are only counted once.
        
        FIXME: Figure out if this can be vectorized somehow instead of using 
        for loops
        '''
        #start = time.time()
        n = matrix.shape[0]
        m = self.lambdifyRegions(matrix)        
        area = 0
        for c in self.domain: #for every point in our domain
            for i in range(n):
                if m[i][0](c) <= m[i][1](c): #if the point is in this region
                    area += 1
                    break #don't need to check the rest of the regions
        #end = time.time()
        #print("time for estimateArea: " + str(end - start))
        return area/((self.numSamples/self.realRange)*(self.numSamples/self.imagRange))
    
    
    def getDomain(self, matrix):
        '''
        getDomain - Automatically calculates the domain for Gershgorwin Regions
        and estimating area. Also stores members realRange and imagRange later
        for scaling in estimateArea()
        @paramter matrix
        @paramter density - an integer specifying the number of sample points
        per unit length
        '''
        #start = time.time()
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
        
        
        maxReal = math.ceil(maxDiagReal + maxRadius)
        minReal = math.floor(minDiagReal - maxRadius)
        maxImag = math.ceil(maxDiagImag + maxRadius)
        minImag = math.floor(minDiagImag - maxRadius)
        
        
        #stored for scaling purposes in estimateArea()
        self.realRange = (maxReal - minReal)
        self.imagRange = (maxImag - minImag)           
            
                         
        domain = np.linspace(minReal, maxReal - 1, self.numSamples)
        domain = np.array([r + 1j*np.linspace(minImag, maxImag - 1, self.numSamples) for r in domain]).flatten()
        #end = time.time()
        #print("time for getDomain: " + str(end - start))
        return domain
    
    
    def getPossibleSets(self, size=4):
        ''' returns a list of lists containing the correct indexes for 
        isospectral reductions. Assuming we're removing 1 index.
        '''
        #start = time.time()
        U = np.arange(0, size)
        possibleSets = []
        for i in range(size):
            possibleSets.append(np.delete(U, i))
            
        #end = time.time()
        #print("time for getPossibleSets: " + str(end - start))
        #print("Possible rows: " + str(possibleRows))
        return possibleSets
    
    
    def getAreas(self):
        '''getAreas - returns a list containing the areas of each reduction
        from a single row. 
        For a 4x4 matrix:
            areas[0] - the area of region obtained by removing row 
        '''
        
        #start = time.time()
        n = self.M.shape[0]
        areas = []
        
        #get original area
        areas.append(self.estimateArea(self.M))
        
        #get areas of possible RSMs
        possibleRows = self.getPossibleSets(size=n)
        for S in possibleRows:
            RSM = self.isospectralReduction(self.M, S)
            areas.append(self.estimateArea(RSM))
        #end = time.time()
        #print("time for getAreas: "  + str(end - start))
        return areas            
    
    
    def syMatrixToNumpyArray(self, syMatrix):
        #start = time.time()
        g = sy.lambdify((), syMatrix, modules="numpy", printer=ImaginaryPrinter)
        #end = time.time()
        #print("time for syMatrixToNumpyArray: " + str(end - start))
        return g()
    
    
    def getIndexSet(self, row, size = 4):
        '''Given a row to remove, getIndexSet returns the corresponding list.
        E.g. if we want to remove row 2 from a 4x4 matrix, we'll return the 
        list [0, 1, 3]. Recall python indwexing starts at 0.
        '''
        U = np.arange(0, size)
        return np.delete(U, row)
    
    
    def getRadii(self, Matrix):
        ''' calculate the radii of Matrix using the sum of the 
        absolute values of the off-diagonal matrix.
        '''
        
        n = Matrix.shape[0]
        radii = []
        for i in range(n):            
            radius = sum([abs(Matrix.row(i)[j]) for j in range(n) if i != j])
            radii.append(radius)
        return radii   
    
