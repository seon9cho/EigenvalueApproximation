# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:06:47 2018

@author: Spencer
"""

"""
Code for the Gershgorwin problem 
"""

import numpy as np
import sympy as sy
import scipy.linalg as la
import matplotlib.pyplot as plt
sy.init_printing()



class Gershgorwin():
    
    def __init__(self, matrix, rowList):
        '''
        Constructor
        
        @parameter matrix - a python list of the matrix we want to use        
        @param rowList - a list specifying the rows of the matrix M to keep. 
        E.g. to keep rows 2, 3, and 4, we would input S = [1, 2, 3]. Recall 
        that python indexing starts at zero.
        
        @member x - a sympy symbol we use in place of lambda
        @member M - the original matrix we're reducing. Stored as a sympy matrix
        @member Mspectrum - spectrum of Matrix M.
        @member S - the rows of M we're keeping in the reduction (indexed at 0)
        @member RSM - the isospectral reduction of M while keeping the rows in S
        @member RSMSpectrum - spectrum of RSM
        @member MGregions - A list of pairs of the Gershgorwin regions of 
        matrix M. Feed this to the plotGRegion function to plot them.
        @member RMSGRegions - A list of pairs of Gershgorwin regions of RSM.
        Feed this to the plotRegion function to plot them. 
        '''
        
        self.x = sy.Symbol('x') 
        self.M = sy.Matrix(matrix)
        self.MSpectrum = la.eigvals(matrix) 
        self.S = rowList
        self.RSM = self.isospectralReduction(rowList) #isospectral reduction of 
        #self.RSMSpectrum = self.RSMSpectrum() #FIXME how to get these into np array
        self.MGRegions = self.gershgorwinRegions(self.M)
        self.RSMGRegions = self.gershgorwinRegions(self.RSM)
        
        self.domain = self.getDomain(matrix, 200)

    
        
    def isospectralReduction(self, rowList):
        '''isospectral Reduction - isospectral reduction of matrix M over the
        over the index set rowList  
        
        @param rowList - a list specifying the rows of the matrix M to keep. E.g. to 
        keep rows 2, 3, and 4, we would input S = [1, 2, 3]. Recall that python
        indexing starts at zero.
        '''
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
        self.RSM = RSM
        return RSM
    
    
    def RSMSpectrum(self):
        character_poly = (self.RSM - self.x * sy.eye(self.RSM.shape[0])).det()
        expr = sy.solve(character_poly, self.x)
        return expr
    
    
    def gershgorwinRegions(self, matrix):
        '''Obtain gershgorwinRegions of matrix
        @ param matrix - the matrix we're geting the gershgorwinRegions for
        '''
    
        n = matrix.shape[0]
        g = []
        for i in range(n):
            lhs = abs(self.x - matrix[i, i])
            rhs = sum([abs(matrix.row(i)[j]) for j in range(n) if i != j])
            expr =  lhs <= rhs
            g.append(expr)
        return g   
               
    
    def plotGRegions(self, matrix):
        '''credit: Seong-Eun
        '''
        n = matrix.shape[0]
        
        #Get a list of tuples of callable functions for the regions
        m = self.lambdifyRegions(matrix)
            
            
        g_range = [[] for i in range(n)]
        for c in self.domain:
            for r, f in zip(g_range, m):
                if f[0](c) <= f[1](c):
                    r.append(c)
        plt.figure(figsize=(10,10))
        for r in g_range:
            plt.scatter(np.real(r), np.imag(r), alpha = 0.1, s = 30)
        plt.scatter(np.real(self.MSpectrum), np.imag(self.MSpectrum), s=50, c='r', label='eigenvalues')
        if matrix == self.M:
            plt.title("M")
        else:
            plt.title("S = " + str(self.S))
        plt.xlabel("real")
        plt.ylabel("imaginary")
        plt.legend()
        plt.axes().set_aspect('equal')
#==============================================================================
#         plt.ylim((-5,5))
#         plt.xlim((-5,5))
#==============================================================================
        plt.show()

        
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
        m = self.lambdifyRegions(matrix)
        area = 0
        for c in self.domain: #for every point in our domain
            for i in range(n):
                if m[i][0](c) <= m[i][1](c): #if the point is in this region
                    area += 1
                    break #don't need to check the rest of the regions          
        return area

        
                    
    def getDomain(self, matrix, numSamples):
        '''
        getDomain - Automatically calculates the domain for Gershgorwin Regions
        and estimating area. 
        @paramter matrix
        @an integer specifying the number of samples to use in each axis
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
        
        #numRealSamples = (maxReal - minReal)*density
        #numImagSamples = (maxImag - minImag)*density 
                         
        numRealSamples = numSamples
        numImagSamples = numSamples                 
                         
        
        domain = np.linspace(minReal, maxReal, numRealSamples)
        domain = np.array([r + 1j*np.linspace(minImag, maxImag, numImagSamples) for r in domain]).flatten()
        self.domain = domain
        return domain


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
    
    
    def getIndexSet(self, row, size = 4):
        '''Given a row to remove, getIndexSet returns the corresponding list.
        E.g. if we want to remove row 2 from a 4x4 matrix, we'll return the 
        list [0, 1, 3]. Recall python indwexing starts at 0.
        '''
        U = np.arange(0, size)
        return np.delete(U, row)






