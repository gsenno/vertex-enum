'''
Created on 27 dic. 2018

@author: gsenno
'''
import numpy as np
import cdd as cdd
from itertools import product,filterfalse
from _functools import reduce
from ncpol2sdpa.sdp_relaxation import imap

class BellPolytope:
    
    outputsAlice = []
    outputsBob = []
    
    # outputsAlice is a list of integers such that outputsAlice[i]==#outputs for Alice's input i (idem outputsBob).
    # the lists do not need to be of the same length.
    def __init__(self,outputsAlice,outputsBob):
        self.outputsAlice = outputsAlice
        self.outputsBob = outputsBob
        
    def _numberOfInputsAlice(self):
        return len(self.outputsAlice)
    
    def _numberOfInputsBob(self):
        return len(self.outputsBob)
    
    def _strategiesGenerator(self,outputs):
        yield from product(*[range(0,numberOfOutputs) for numberOfOutputs in outputs])
    
    def _strategyToDistribution(self,stgAlice, stgBob):
        distribution = []
        for x in range (0,len(self.outputsAlice)):
            for y in range (0,len(self.outputsBob)):
                for a in range (0,self.outputsAlice[x]):
                    for b in range (0,self.outputsBob[y]):
                        if (a==stgAlice[x])&(b==stgBob[y]):
                            distribution.append(1)
                        else:
                            distribution.append(0)
        return distribution
    
    def getGeneratorForVertices(self):
        return (self._strategyToDistribution(stgAlice, stgBob)
                 for stgAlice in self._strategiesGenerator(self.outputsAlice)
                 for stgBob in self._strategiesGenerator(self.outputsBob)) 
    
    def getListOfVertices(self):
        return list(self.getGeneratorForVertices())
    
    
