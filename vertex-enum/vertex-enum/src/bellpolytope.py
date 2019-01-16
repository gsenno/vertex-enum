'''
Created on 27 dic. 2018

@author: gsenno
'''
import numpy as np
import cdd as cdd
from itertools import product

class BellPolytope:
    parties = 2
    N = 0
    K = 0
    vertices = []
    inequalities = []
    
    def __init__(self,inputsPerParty,outputsPerParty):
        self.K = outputsPerParty
        self.N = inputsPerParty
    
    def getVertices(self):  # @DontTrace
        if self.vertices==[]:
            self.vertices=self.__generateVertices(self.parties,self.N,self.K)
        return self.vertices
    
    def numberOfVertices(self):
        return len(self.getVertices())
    
    def getInequalities(self):  # @DontTrace
        if self.inequalities==[]:
            self.inequalities=self.__generateInequalities(self.getVertices())
        return self.inequalities
    
    def getInefficiencyResistantInequalities(self):
        origen = np.zeros((self.K*self.N)**self.parties)
        verticesLocalIncomplete = list(self.getVertices())+[origen];
    
        inequalities = self.__generateInequalities(verticesLocalIncomplete)
    
        unNormalizedIneffResistIneq = map(lambda ineq: self.__extendInequalityToDetecLoopholeSetting(ineq),inequalities);
    
        ineffResistInequalities = map(lambda ineq: self.__makeIneqBoundedOnAbortDist(ineq),unNormalizedIneffResistIneq)
        
        return ineffResistInequalities
    
    def __generateVertices(self,parties,inputsPerParty,outputsPerParty):
        D=np.zeros((outputsPerParty**inputsPerParty,inputsPerParty), dtype=int)
        for _ in range(outputsPerParty**inputsPerParty):
            D[_][:]=np.array(np.unravel_index(_,(outputsPerParty,)*inputsPerParty))
        vertices=np.zeros(((outputsPerParty**(inputsPerParty*parties),)+(inputsPerParty,)*parties+(outputsPerParty,)*parties))
        c=0
        for _ in product(range(outputsPerParty**inputsPerParty), repeat=parties):
            for x in product(range(inputsPerParty), repeat=parties):
                vertices[(c,)+x+tuple([D[_[i]][x[i]] for i in range(parties)])]=1
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
        inputs=self.N**self.parties
        outputs=self.K**self.parties
        ineffResistInequality = inequality[0:1]
        for nInputPair in range(inputs):
            oldCoeffForInputs=functional[nInputPair*outputs:(nInputPair+1)*outputs]
            newCoeffForInputs=np.zeros((self.K+1)**self.parties)
            for nOutput in range(self.K):
                newCoeffForInputs[(self.K+1)*nOutput:(self.K+1)*nOutput+self.K] = oldCoeffForInputs[nOutput*self.K:nOutput*self.K+self.K] 
            ineffResistInequality=np.append(ineffResistInequality,newCoeffForInputs)
        
        return ineffResistInequality
    
    def __makeIneqBoundedOnAbortDist(self,ineq):
        abortingDists = self.__generateVertices(self.parties,self.N,self.K+1)
        bound = max(map(lambda vector : np.dot(ineq[1:],vector),abortingDists))
        ineffResistIneq = np.concatenate(([0],np.array(ineq[1:]))) if bound==0 else np.concatenate(([1],(1/bound)*np.array(ineq[1:]))) 
        return ineffResistIneq
