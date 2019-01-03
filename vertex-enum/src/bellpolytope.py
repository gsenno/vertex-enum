'''
Created on 27 dic. 2018

@author: gsenno
'''
import numpy as np
import cdd as cdd
from itertools import product, permutations, combinations, combinations_with_replacement

class BellPolytope:
    K = 0
    N = 0
    vertices = []
    inequalities = []
    
    def __init__(self,inputsPerParty,outputsPerParty):
        self.K = outputsPerParty
        self.N = inputsPerParty
    
    def getVertices(self):  # @DontTrace
        if self.vertices==[]:
            self.vertices=self.__generateVertices(self.K,self.N)
        return self.vertices
    
    def numberOfVertices(self):
        return len(self.getVertices())
    
    def getInequalities(self):  # @DontTrace
        if self.inequalities==[]:
            self.inequalities=self.__generateInequalities(self.getVertices())
        return self.inequalities
    
    def getInefficiencyResistantInequalities(self):
        origen = np.zeros((self.K*self.N)**2)
        verticesLocalIncomplete = list(self.getVertices())+[origen];
    
        inequalities = self.__generateInequalities(verticesLocalIncomplete)
    
        unNormalizedIneffResistIneq = map(lambda ineq: self.__extendInequalityToDetecLoopholeSetting(ineq),inequalities);
    
        ineffResistInequalities = map(lambda ineq: self.__makeIneqBoundedOnAbortDist(ineq),unNormalizedIneffResistIneq)
        
        return ineffResistInequalities
    
    def __generateVertices(self,K,N):
        D=np.zeros((K**N,N), dtype=int)
        for _ in range(K**N):
            D[_][:]=np.array(np.unravel_index(_,(K,)*N))
        vertices=np.zeros(((K**(N*2),)+(N,)*2+(K,)*2))
        c=0
        for _ in product(range(K**N), repeat=2):
            for x in product(range(N), repeat=2):
                vertices[(c,)+x+tuple([D[_[i]][x[i]] for i in range(2)])]=1
            c+=1
        shape=np.prod(np.delete(vertices.shape[:],0,0))
        return vertices.reshape(vertices.shape[0],shape)
    
    def __generateInequalities(self,vertices):
        cddPolytope = self.__generateCddPolytope(vertices)
        ext = cddPolytope.get_inequalities()
        inequalities = map(lambda x:list(x), list(ext.__getitem__(slice(ext.row_size))))
        return inequalities
    
    def __generateCddPolytope(self,vertices):
        vRep = self.__buildVRepresentation(vertices)
        mat = cdd.Matrix(vRep, number_type='fraction')
        mat.rep_type = cdd.RepType.GENERATOR
        poly = cdd.Polyhedron(mat)
        return poly
    
    def __buildVRepresentation(self,vertices):
        vRep = np.zeros((len(vertices), 1+len(vertices[0])))
        for i in range(len(vertices)):
            vRep[i][0] = 1
            vRep[i][1:] = vertices[i]
        
        return vRep
    
    def __extendInequalityToDetecLoopholeSetting(self,inequality):
        functional=inequality[1:]
        inputs=self.N**2
        outputs=self.K**2
        ineffResistInequality = inequality[0:1]
        for nInputPair in range(inputs):
            oldCoeffForInputs=functional[nInputPair*outputs:(nInputPair+1)*outputs]
            newCoeffForInputs=np.zeros((self.K+1)**2)
            for nOutput in range(self.K):
                newCoeffForInputs[(self.K+1)*nOutput:(self.K+1)*nOutput+self.K] = oldCoeffForInputs[nOutput*self.K:nOutput*self.K+self.K] 
            ineffResistInequality.extend(newCoeffForInputs)
        
        return ineffResistInequality
    
    def __makeIneqBoundedOnAbortDist(self,ineq):
        abortingDists = self.__generateVertices(self.K+1, self.N)
        bound = max(map(lambda vector : np.dot(ineq[1:],vector),abortingDists))
        ineffResistIneq = np.concatenate(([0],np.array(ineq[1:]))) if bound==0 else np.concatenate(([1],(1/bound)*np.array(ineq[1:]))) 
        return ineffResistIneq
