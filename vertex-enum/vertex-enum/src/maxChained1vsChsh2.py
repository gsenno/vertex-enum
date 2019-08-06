
'''
Created on Feb 12, 2019

@author: rravell
'''
from mosek.fusion import *
import numpy as np
import cvxopt as cvx
import itertools as it
from ncpol2sdpa.sdp_relaxation import imap
from linopttools import *
import qutip as qt
from itertools import product, islice
import picos as pic
from bellpolytopewithonewaycomm import BellPolytopeWithOneWayCommunication

def CHAINED(n,A,B):
    result = 0
    for i in range(0,n-1):
        result+=pic.kron(A[i],B[i])+pic.kron(A[i+1],B[i])
    result+=pic.kron(A[n-1],B[n-1])-pic.kron(A[0],B[n-1])
    return result

def CHSH(A,B):
    #chsh bell op give observables
    return pic.kron(A[0],B[0])+pic.kron(A[0],B[1])+pic.kron(A[1],B[0]) \
            -pic.kron(A[1],B[1])


def chainedBellValue(n,p):
    result = 0
    for i in range(n-1):
        for a,b in product(range(2),repeat=2):
            result+=(-1)**(a+b)*(p[(b+2*a)+4*(i+n*i)]+p[(b+2*a)+4*(i+n*(i+1))])
    for a,b in product(range(2),repeat=2):
        result+=(-1)**(a+b)*(p[(b+2*a)+4*(n-1+n*(n-1))]-p[(b+2*a)+4*(n-1+n*(0))])
    return result

if __name__ == '__main__':
    
    n=3
    outputsAlice = [4,4,4,4,4,4]
    outputsBob = [2 for _ in range(0,n)]
               
    alpha=5.1
     
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
     
    prob=pic.Problem() 
     
    A={}
    for x1,x2 in product(range(n),range(2)):
        for a1,a2 in product(range(2),repeat=2):
            A[x1,x2,a1,a2]=prob.add_variable('A_{0}{1}{2}{3}'.format(x1,x2,a1,a2),
                                             (2,2),'hermitian')
            prob.add_constraint(A[x1,x2,a1,a2]>>0)
     
    for x1,x2 in product(range(n),range(2)):
        prob.add_constraint(sum(A[x1,x2,a1,a2] for a1 in range(2) 
                                               for a2 in range(2))==np.eye(2))
    
    phi = lambda i : (i-1)*np.pi/n
    phiprime= lambda i : (2*i-1)*np.pi/(2*n) 
    bobEffects=[projectorsForQubitObservable
                    (createQubitObservable([np.sin(phiprime(i)),0,np.cos(phiprime(i))]))
                     for i in range(1,n+1)]
    bobEffects=[list(map(lambda qutipProj : qutipProj.get_data().toarray(),obs)) for obs in bobEffects]
     
     
    B={}
    for i in range(n):
        for j in (0,1):
            B[i,j]=pic.new_param('B_'+str(i)+str(j),bobEffects[i][j])
     
    rho=pic.new_param('rho',np.outer([1,0,0,1],[1,0,0,1])/2)
     
     
    A1=[1/2*sum(A[x1,x2,a1,a2]*(-1)**a1 for a1 in range(2) 
                                        for a2 in range(2) 
                                        for x2 in range(2))
                                                    for x1 in range(n)]
     
    A2=[1/n*sum(A[x1,x2,a1,a2]*(-1)**a2 for a1 in range(2) 
                                        for a2 in range(2) 
                                        for x1 in range(n))
                                                    for x2 in range(2)]
         
    B1=[sum(B[y,b]*(-1)**b for b in range(2)) for y in range(n)]
     
     
    prob.add_constraint(pic.trace(CHAINED(n,A1,B1)*rho)==alpha)
     
    prob.set_objective('max',
                       pic.trace(CHSH(A2,B1)*rho))
     
    prob.solve()
#     
    dist=[pic.trace(pic.kron(A[x1,x2,a1,a2],B[y,b])*rho).get_value().real
          for x1,x2,y,a1,a2,b in product(range(n),range(2),range(n),range(2),range(2),range(2))]
 
    vertices=BellPolytopeWithOneWayCommunication(outputsAlice,outputsBob).getListOfVertices()
    
    with Model("lo1") as M:

        # Create variable 'x' of length 4
        x = M.variable("x", len(vertices), Domain.greaterThan(0.0))
        

        # Create constraints
        for prob in range(len(vertices[0])):
            M.constraint('p('+str(prob)+')',Expr.dot(x,list(map(lambda ver : ver[prob],vertices))),Domain.equalsTo(dist[prob]))
            
        M.constraint('norm',Expr.dot(x,np.ones((len(vertices), 1))),Domain.equalsTo(1))
        
        
        # Set the objective function to (c^t * x)
        M.objective("obj", ObjectiveSense.Maximize, 1)

        # Solve the problem
        M.solve()

        # Get the solution values
        print(M.getProblemStatus(SolutionType.Basic))
        print(x.dual)
