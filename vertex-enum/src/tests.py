import numpy as np
import unittest
from bellpolytope import BellPolytope

class TestBellPolytope(unittest.TestCase):

#     def testIneqsAreIneffResist(self):    
#         inputs = 3
#         outputs = 2
#         poly = BellPolytope(3,2)    
#         ineffResistInequalities = poly.getInefficiencyResistantInequalities()
#         
#         K = outputs+1
#         totalOutputs = K**2
#         
#         N=inputs
#         totalInputs = N**2
#         
#         nZeroConstraints = 2*totalInputs*(totalOutputs-(K-1)**2)
#         nLowerBoundOnDeterministicConstraints = (K**(2*N))
#         nUpperBoundOnDeterministicConstraints = (K**(2*N))
#         nConstraints =  nZeroConstraints + nLowerBoundOnDeterministicConstraints + nUpperBoundOnDeterministicConstraints
#         nVariables = (K*N)**2
#         
#         matrixHRepresentation=np.zeros((nConstraints,1+nVariables))
#         
#         #IMPOSING THE COMPONENTS 0 OF THE FUNCTIONAL
#         rowNum = 0
#         for j in range (1, 1+K*(totalInputs)):
#             matrixHRepresentation[j-1, K*j]=1
#             rowNum+=1
#               
#         for j in range (0, totalInputs):
#             matrixHRepresentation[rowNum + j, (K-1)*K+1 + totalOutputs*j]= 1
#             matrixHRepresentation[rowNum + 1 + j, (K-1)*K+2 + totalOutputs*j] = 1
#             rowNum+=2 
#          
#         #VALUES OF THE INEQUALITIES
#         matrixHRepresentation[rowNum:rowNum+nLowerBoundOnDeterministicConstraints,0]=63
#         matrixHRepresentation[rowNum+nLowerBoundOnDeterministicConstraints:rowNum+nLowerBoundOnDeterministicConstraints+nUpperBoundOnDeterministicConstraints,0]=1 
#         
#         i=0
#         
#         poly = BellPolytope(K,N)
#         for vector in poly.getVertices():
#             matrixHRepresentation[rowNum+i,1:]=vector
#             matrixHRepresentation[rowNum+nLowerBoundOnDeterministicConstraints+i,1:]=-vector
#             i+=1
#         
#         for j in range (1,1+K*(totalInputs)):
#             matrixHRepresentation[rowNum+j,K*j]=-1
#             rowNum+=1
#         
#         for j in range (0, totalInputs):
#             matrixHRepresentation[rowNum + j, (K-1)*K+1 + totalOutputs*j]= -1
#             matrixHRepresentation[rowNum + 1 + j, (K-1)*K+2 + totalOutputs*j] = -1
#             rowNum+=2                         
#                                  
#                                
#         
#         mat=matrixHRepresentation[:,1:]
#         b=matrixHRepresentation[:,0]
#         
#         
#         for ineq in ineffResistInequalities:
#             results = np.greater_equal(np.dot(mat,ineq[1:]),(-1)*b)
#             map(lambda res : self.assertTrue(res, 'ineq '+np.array2string(ineq)+' is not ineff resist'),results)
#         
    
    def testNumberOfVerticesOfBellPolytope23Is81(self):
        self.assertEqual(BellPolytope(2,3).numberOfVertices(),3**4,'')
    
    def testNumberOfVerticesOfBellPolytope33Is81(self):
        self.assertEqual(BellPolytope(3,3).numberOfVertices(),3**6,'')
        
    def testAliceAndBobAlways1IsVertexOfBellPolytope33(self):
        pass
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                            
                            
                            
                        