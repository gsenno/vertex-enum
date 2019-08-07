
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
import picos as pic

def In_ConvexHull(vertices,q):
    # Tests if the point q is inside the convex Hull of the points D
    # q should be a np multidim array
    # D should be any np multidim array with first index labelling the points of convex set
    # output is list containing the solver status and the solution 
    #reshape so we have vectors
    
    N=len(vertices)
    D=pic.new_param('D',vertices)
    #define problem
    prob=pic.Problem()
    #cerate prob vector
    p=prob.add_variable('p',N)
    #add desired point
    q=np.reshape(q,[1,-1])
    q=pic.new_param('q',q)
    #feasibilitiy test
    prob.set_objective('max',0*p[0])
    
    #constraints: positivity, normalisation, correct vector
    prob.add_constraint(p>=0)
    prob.add_constraint([[1 for __ in range(N)]]*p==1)
#    prob.add_constraint(pic.sum([p[i] for i in range(N)])==1.0)
    prob.add_constraint(p.T*D==q)

    prob.solve(verbose=1)
#    print(prob)
    
    #get optimal variables and reshape
    pvalues=np.array(p.value)
    
    return [prob.status, pvalues]

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
     
    print(In_ConvexHull(vertices, vertices[0]))
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
