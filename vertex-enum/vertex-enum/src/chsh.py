
'''
Created on Feb 12, 2019

@author: rravell
'''
from scipy.spatial import ConvexHull
from mosek.fusion import *
from bellpolytope import BellPolytope
import cdd as cdd
import numpy as np
import itertools as it
from _functools import reduce
from ncpol2sdpa.sdp_relaxation import imap
from linopttools import *
import qutip as qt

if __name__ == '__main__':
    
    outputsAlice = [2,2]
    outputsBob = [2,2]
    
    aliceUnBlochVectors = [[1,0,0],[0,1,0]]
    aliceObservables = list(map(lambda bloch : createQubitObservable(bloch),aliceUnBlochVectors))
    aliceEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),aliceObservables))
    
    bobUnBlochVectors = [[-1,-1,0],[-1,1,0]]
    bobObservables=list(map(lambda bloch : createQubitObservable(bloch),bobUnBlochVectors))
    bobEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),bobObservables))
    
    psi=createMaxEntState(2)
    
    dist=computeDistributionFromStateAndEffects(psi,aliceEffects,bobEffects)
 
    vertices=BellPolytope(outputsAlice,outputsBob).getListOfVertices()
    
    ConvexHull(vertices).
#     
#     with Model("lo1") as M:
#         M.setSolverParam('logInfeasAna',100)
#         # Create variable 'x' of length 4
#         x = M.variable("x", len(vertices), Domain.greaterThan(0.0))
# 
#         # Create constraints
#         constraints=[]
#         for prob in range(len(vertices[0])):
#             constraints.append(M.constraint('p('+str(prob)+')',Expr.dot(x,list(map(lambda ver : ver[prob],vertices))),Domain.equalsTo(dist[prob])))
#             
#         constraints.append(M.constraint('norm',Expr.dot(x,np.ones((len(vertices), 1))),Domain.equalsTo(1)))
#         
#         
#         # Set the objective function to (c^t * x)
#         M.objective("obj", ObjectiveSense.Minimize, 1)
# 
#         # Solve the problem
#         report = M.solve()
# 
#         # Get the solution values
#         print(M.getProblemStatus(SolutionType.Basic))
#         
