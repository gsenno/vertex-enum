'''
Created on 5 ago. 2019

@author: gsenno
'''
from bellpolytope import BellPolytope
from itertools import product
from _functools import reduce


#Polytope for bipartite strategies using shared randomness and communication from Alice to Bob.
class BellPolytopeWithOneWayCommunication(BellPolytope):
    '''
    classdocs
    '''

    #see bellpolytope.py for description of params
    def __init__(self, outputsAlice, outputsBob):
        BellPolytope.__init__(self, outputsAlice, outputsBob)

    def getGeneratorForVertices(self):
        #local vertices
        yield from BellPolytope.getGeneratorForVertices(self)
        
        #distributions with nontrivial use of the communication channel
        communicationStrgs=[format(i,'0'+str(self._numberOfInputsAlice())+'b') 
                            for i in range(1,2**(self._numberOfInputsAlice()-1))]
        strgsAlice = [[(stgAlice[i],int(comm[i])) for i in range(0,len(stgAlice))] 
                      for stgAlice in self._strategiesGenerator(self.outputsAlice) for comm in communicationStrgs]
        
        strgsBob = [stgBob for stgBob in 
                        self._strategiesGenerator(reduce(lambda acum,elem : acum+[elem,elem],self.outputsBob,[]))
                    if stgBob[0::2]!=stgBob[1::2]]
        
        yield from ([int(a==stgAlice[x][0])&(b==stgBob[2*y+stgAlice[x][1]])
                     for x in range(self._numberOfInputsAlice())
                     for y in range(self._numberOfInputsBob())
                     for a in range(self.outputsAlice[x])
                     for b in range(self.outputsBob[y])] 
               for stgAlice in strgsAlice for stgBob in strgsBob)
